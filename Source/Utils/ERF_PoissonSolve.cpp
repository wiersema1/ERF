#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

/**
 * Define the domain boundary conditions for the (optional) Poisson solve
 * if we want to enforce incompressibility of the initial conditions
 */

using BCType = LinOpBCType;

Array<LinOpBCType,AMREX_SPACEDIM>
ERF::get_projection_bc (Orientation::Side side) const noexcept
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc_type = domain_bc_type[Orientation(dir,side)];
            if (bc_type == "Outflow") {
                r[dir] = LinOpBCType::Dirichlet;
            } else
            {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    return r;
}
bool ERF::projection_has_dirichlet (Array<LinOpBCType,AMREX_SPACEDIM> bcs) const
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (bcs[dir] == LinOpBCType::Dirichlet) return true;
    }
    return false;
}


/**
 * Project the single-level velocity field to enforce incompressibility
 * Note that the level may or may not be level 0.
 */
void ERF::project_velocities (int lev, Real l_dt, Vector<MultiFab>& mom_mf, MultiFab& /*Omega*/, MultiFab& pmf)
{
    BL_PROFILE("ERF::project_velocities()");

    AMREX_ALWAYS_ASSERT(!solverChoice.use_terrain);

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    LPInfo info;
    // Allow a hidden direction if the domain is one cell wide in any lateral direction
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
    }

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(mom_mf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(mom_mf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    MLPoisson mlpoisson(geom_tmp, ba_tmp, dm_tmp, info);

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);
    bool need_adjust_rhs = (projection_has_dirichlet(bclo) || projection_has_dirichlet(bchi)) ? false : true;
    mlpoisson.setDomainBC(bclo, bchi);

    if (lev > 0) {
        mlpoisson.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
    mlpoisson.setLevelBC(0, nullptr);

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;

    rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
    phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1);

    fluxes.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
    }

    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;
    rho0_u_const[0] = &mom_mf[IntVars::xmom];
    rho0_u_const[1] = &mom_mf[IntVars::ymom];
    rho0_u_const[2] = &mom_mf[IntVars::zmom];

    computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);

    if (mg_verbose > 0) {
        Print() << "Max norm of divergence before at level " << lev << " : " << rhs[0].norm0() << std::endl;
    }

    // If all Neumann BCs, adjust RHS to make sure we can converge
    if (need_adjust_rhs)
    {
        Real offset = volWgtSumMF(lev, rhs[0], 0, *mapfac_m[lev], false, false);
        if (mg_verbose > 1) {
            Print() << "Poisson solvability offset = " << offset << std::endl;
        }
        rhs[0].plus(-offset, 0, 1);
    }

    Real start_step = static_cast<Real>(ParallelDescriptor::second());

#ifdef ERF_USE_FFT
    if (use_fft) {
        AMREX_ALWAYS_ASSERT(lev == 0);
        if (!m_poisson) {
            m_poisson = std::make_unique<FFT::PoissonHybrid<MultiFab>>(Geom(0));
        }
        m_poisson->solve(phi[lev], rhs[lev]);

        phi[lev].FillBoundary(geom[lev].periodicity());

        auto dxInv = geom[lev].InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const&  p_arr  = phi[lev].array(mfi);

        Box const& xbx = mfi.nodaltilebox(0);
        const Real dx_inv = dxInv[0];
        Array4<Real> const& fx_arr  = fluxes[lev][0].array(mfi);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i-1,j,k)) * dx_inv;
        });

        Box const& ybx = mfi.nodaltilebox(1);
        const Real dy_inv = dxInv[1];
        Array4<Real> const& fy_arr  = fluxes[lev][1].array(mfi);
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j-1,k)) * dy_inv;
        });

        auto const dom_lo = lbound(geom[lev].Domain());
        auto const dom_hi = ubound(geom[lev].Domain());

        Box const& zbx = mfi.nodaltilebox(2);
        const Real dz_inv = dxInv[2];
        Array4<Real> const& fz_arr  = fluxes[lev][2].array(mfi);
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k == dom_lo.z || k == dom_hi.z+1) {
                fz_arr(i,j,k) = 0.0;
            } else {
                fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j,k-1)) * dz_inv;
            }
        });
    } // mfi
     } else
#endif
    {
        // Initialize phi to 0
        phi[0].setVal(0.0);

        MLMG mlmg(mlpoisson);
        int max_iter = 100;
        mlmg.setMaxIter(max_iter);

        mlmg.setVerbose(mg_verbose);
        mlmg.setBottomVerbose(0);

        mlmg.solve(GetVecOfPtrs(phi),
                   GetVecOfConstPtrs(rhs),
                   solverChoice.poisson_reltol,
                   solverChoice.poisson_abstol);
        mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));
    }

    // Subtract dt grad(phi) from the momenta
    MultiFab::Add(mom_mf[IntVars::xmom],fluxes[0][0],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::ymom],fluxes[0][1],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::zmom],fluxes[0][2],0,0,1,0);

    // Update pressure variable with phi -- note that phi is dt * change in pressure
    MultiFab::Saxpy(pmf, 1.0/l_dt, phi[0],0,0,1,0);
    pmf.FillBoundary(geom[lev].periodicity());

    Real end_step = static_cast<Real>(ParallelDescriptor::second());
    if (mg_verbose > 0) {
        amrex::Print() << "Time in solve " << end_step - start_step << std::endl;
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mom_mf[Vars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& pp_arr  = pmf.array(mfi);
        Box const& bx    = mfi.tilebox();
        auto const bx_lo = lbound(bx);
        auto const bx_hi = ubound(bx);
        if (bx_lo.x == dom_lo.x) {
            if (bclo[0] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,0,dom_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i-1,j,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,0,dom_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i-1,j,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_lo.y == dom_lo.y) {
            if (bclo[1] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_lo.z == dom_lo.z) {
            ParallelFor(makeSlab(bx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k-1) = pp_arr(i,j,k);
            });
        }
        if (bx_hi.x == dom_hi.x) {
            if (bchi[0] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,0,dom_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i+1,j,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,0,dom_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i+1,j,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_hi.y == dom_hi.y) {
            if (bchi[1] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_hi.z == dom_hi.z) {
            ParallelFor(makeSlab(bx,2,dom_hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k+1) = pp_arr(i,j,k);
            });
        }
    } // mfi

    // Now overwrite with periodic fill outside domain and fine-fine fill inside
    pmf.FillBoundary(geom[lev].periodicity());

    //
    // BELOW IS SIMPLY VERIFYING THE DIVERGENCE AFTER THE SOLVE
    //
    if (mg_verbose > 0) {
        computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);
        Print() << "Max norm of divergence after solve at level " << lev << " : " << rhs[0].norm0() << std::endl;
    }
}
