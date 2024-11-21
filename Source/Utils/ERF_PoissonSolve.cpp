#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Project the single-level velocity field to enforce incompressibility
 * Note that the level may or may not be level 0.
 */
void ERF::project_velocities (int lev, Real l_dt, Vector<MultiFab>& mom_mf, MultiFab& pmf)
{
    BL_PROFILE("ERF::project_velocities()");

    bool l_use_terrain = SolverChoice::terrain_type != TerrainType::None;
    bool use_gmres     = (l_use_terrain && !SolverChoice::terrain_is_flat);

#ifndef ERF_USE_EB
    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());
#endif

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(mom_mf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(mom_mf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;

#ifdef ERF_USE_EB
    rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0, MFInfo(), Factory(lev));
    phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1, MFInfo(), Factory(lev));
#else
    rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
    phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1);
#endif

    fluxes.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
#ifdef ERF_USE_EB
        fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0, MFInfo(), Factory(lev));
#else
        fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
#endif
    }

    auto dxInv = geom[lev].InvCellSizeArray();

    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;
    rho0_u_const[0] = &mom_mf[IntVars::xmom];
    rho0_u_const[1] = &mom_mf[IntVars::ymom];
    rho0_u_const[2] = &mom_mf[IntVars::zmom];

    Real abstol = solverChoice.poisson_abstol;

    // ****************************************************************************
    // Compute divergence which will form RHS
    // Note that we replace "rho0w" with the contravariant momentum, Omega
    // ****************************************************************************
#ifdef ERF_USE_EB
    bool already_on_centroids = true;
    EB_computeDivergence(rhs[0], rho0_u_const, geom_tmp[0], already_on_centroids);
#else
    if (l_use_terrain)
    {
        for ( MFIter mfi(rhs[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Array4<Real const>& rho0u_arr = mom_mf[IntVars::xmom].const_array(mfi);
            const Array4<Real const>& rho0v_arr = mom_mf[IntVars::ymom].const_array(mfi);
            const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);

            const Array4<Real const>&     z_nd = z_phys_nd[lev]->const_array(mfi);
            //
            // Define Omega from (rho0 W) but store it in the same array
            //
            Box tbz = mfi.nodaltilebox(2);
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (k > dom_lo.z && k <= dom_hi.z) {
                    Real rho0w = rho0w_arr(i,j,k);
                    rho0w_arr(i,j,k) = OmegaFromW(i,j,k,rho0w,rho0u_arr,rho0v_arr,z_nd,dxInv);
                } else {
                    rho0w_arr(i,j,k) = Real(0.0);
                }
            });
        } // mfi
        //
        // Note we compute the divergence after we convert rho0W --> Omega
        //
        for ( MFIter mfi(rhs[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            const Array4<Real const>& rho0u_arr = mom_mf[IntVars::xmom].const_array(mfi);
            const Array4<Real const>& rho0v_arr = mom_mf[IntVars::ymom].const_array(mfi);
            const Array4<Real const>& rho0w_arr = mom_mf[IntVars::zmom].const_array(mfi);
            const Array4<Real      >&  rhs_arr = rhs[0].array(mfi);

            Real* stretched_dz_d_ptr = stretched_dz_d[lev].data();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real dz = stretched_dz_d_ptr[k];
                rhs_arr(i,j,k) =  (rho0u_arr(i+1,j,k) - rho0u_arr(i,j,k)) * dxInv[0]
                                 +(rho0v_arr(i,j+1,k) - rho0v_arr(i,j,k)) * dxInv[1]
                                 +(rho0w_arr(i,j,k+1) - rho0w_arr(i,j,k)) / dz;
            });
        } // mfi

    } else {

        computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);

    }
#endif

    Real rhsnorm = rhs[0].norm0();

    if (mg_verbose > 0) {
        Print() << "Max norm of divergence before at level " << lev << " : " << rhsnorm << std::endl;
    }

    // Initialize phi to 0
    // (It is essential that we do this in order to fill the corners; these are never
    //  used but the Saxpy requires the values to be initialized.)
    phi[0].setVal(0.0);

    if (rhsnorm <= abstol) return;

    Real start_step = static_cast<Real>(ParallelDescriptor::second());

    // ****************************************************************************
    // Choose the solver and solve
    // ****************************************************************************
#ifdef ERF_USE_EB
    solve_with_EB_mlmg(lev, rhs, phi, fluxes);
#else
#ifdef ERF_USE_FFT
    if (use_fft) {

        solve_with_fft(lev, rhs[0], phi[0], fluxes[0]);

    } else
#endif
    if (use_gmres)
    {
        solve_with_gmres(lev, rhs, phi, fluxes);

    } else {

        solve_with_mlmg(lev, rhs, phi, fluxes);
    }
#endif

    // ****************************************************************************
    // Subtract dt grad(phi) from the momenta (rho0u, rho0v, Omega)
    // ****************************************************************************
    MultiFab::Add(mom_mf[IntVars::xmom],fluxes[0][0],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::ymom],fluxes[0][1],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::zmom],fluxes[0][2],0,0,1,0);

    // ****************************************************************************
    // Print time in solve
    // ****************************************************************************
    Real end_step = static_cast<Real>(ParallelDescriptor::second());
    if (mg_verbose > 0) {
        amrex::Print() << "Time in solve " << end_step - start_step << std::endl;
    }

    //
    // This call is only to verify the divergence after the solve
    // It is important we do this before computing the rho0w_arr from Omega back to rho0w
    //
    //
    // ****************************************************************************
    // THIS IS SIMPLY VERIFYING THE DIVERGENCE AFTER THE SOLVE
    // ****************************************************************************
    //
    if (mg_verbose > 0)
    {
        if (l_use_terrain)
        {
            //
            // Note we compute the divergence before we convert Omega back to rho0W
            //
            for ( MFIter mfi(rhs[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                const Array4<Real const>& rho0u_arr = mom_mf[IntVars::xmom].const_array(mfi);
                const Array4<Real const>& rho0v_arr = mom_mf[IntVars::ymom].const_array(mfi);
                const Array4<Real const>& rho0w_arr = mom_mf[IntVars::zmom].const_array(mfi);
                const Array4<Real      >&  rhs_arr = rhs[0].array(mfi);

                Real* stretched_dz_d_ptr = stretched_dz_d[lev].data();
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    Real dz = stretched_dz_d_ptr[k];
                    rhs_arr(i,j,k) =  (rho0u_arr(i+1,j,k) - rho0u_arr(i,j,k)) * dxInv[0]
                                     +(rho0v_arr(i,j+1,k) - rho0v_arr(i,j,k)) * dxInv[1]
                                     +(rho0w_arr(i,j,k+1) - rho0w_arr(i,j,k)) / dz;
                });
            } // mfi

        } else {
            computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);
        }

        amrex::Print() << "Max norm of divergence after solve at level " << lev << " : " << rhs[0].norm0() << std::endl;
    } // mg_verbose

    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab back to holding (rho0w) rather than Omega
    // ****************************************************************************
    //
    if (l_use_terrain)
    {
        for (MFIter mfi(mom_mf[Vars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
             Box tbz = mfi.nodaltilebox(2);
             const Array4<Real      >& rho0u_arr = mom_mf[IntVars::xmom].array(mfi);
             const Array4<Real      >& rho0v_arr = mom_mf[IntVars::ymom].array(mfi);
             const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);
             const Array4<Real const>&      z_nd = z_phys_nd[lev]->const_array(mfi);
             ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                 Real omega = rho0w_arr(i,j,k);
                 rho0w_arr(i,j,k) = WFromOmega(i,j,k,omega,rho0u_arr,rho0v_arr,z_nd,dxInv);
             });
        } // mfi
    }

    // ****************************************************************************
    // Update pressure variable with phi -- note that phi is dt * change in pressure
    // ****************************************************************************
    MultiFab::Saxpy(pmf, 1.0/l_dt, phi[0],0,0,1,1);
    pmf.FillBoundary(geom[lev].periodicity());
}
