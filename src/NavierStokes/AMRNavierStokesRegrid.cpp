/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/somarhub.
 ******************************************************************************/
#include "AMRNavierStokes.H"
#include "AMRNSF_F.H"
#include "StratUtils.H"
#include "MappedFineInterp.H"
#include "computeMappedNorm.H"
#include "MappedAMRPoissonOpFactory.H"
#include "BoxIterator.H"
#include "AMRCCProjector.H"
#include "SetValLevel.H"
#include "ExtrapolationUtils.H"
#include <iomanip>


// -----------------------------------------------------------------------------
// Create a load-balanced DisjointBoxLayout from a collection of boxes
// -----------------------------------------------------------------------------
DisjointBoxLayout AMRNavierStokes::loadBalance (const Vector<Box>& a_grids)
{
    if (s_verbosity >= 6) {
        pout () << "AMRNavierStokes::loadBalance " << m_level
                << "  boxes=" << a_grids.size() << endl;
    }

    Vector<int> proc_map;
    LoadBalance(proc_map, a_grids);

    if (s_verbosity >= 4) {
        pout () << "AMRNavierStokes::loadBalance: processor map: " << endl;
        for (int igrid=0; igrid< a_grids.size(); ++igrid) {
            pout() << igrid << ": " << proc_map[igrid] << "\t" << a_grids[igrid] << endl;
        }
        pout() << endl;
    }

    return DisjointBoxLayout(a_grids, proc_map, m_problem_domain);
}


// -----------------------------------------------------------------------------
// Create tags at initial time.
// The only data available are the initial velocity and scalars --
// basically, whatever was set in the initialData() function.
// -----------------------------------------------------------------------------
void AMRNavierStokes::tagCellsInit (IntVectSet& a_tags)
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::tagCellsInit " << m_level << endl;
    }

    // Chombo has some s_specifyInitialGrids code here. I don't think this is
    // the right place for all that. These tags are merely suggestions for
    // the MeshRefine object, not the end grids themselves.

    // Use same function at initial time as in the rest of the time
    this->tagCells(a_tags);
}


// -----------------------------------------------------------------------------
// Tag cells for regridding
// -----------------------------------------------------------------------------
void AMRNavierStokes::tagCells (IntVectSet& a_tags)
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::tagCells " << m_level << endl;
    }
    CH_TIME("AMRNavierStokes::tagCells");

    const DisjointBoxLayout& grids = newVel().getBoxes();
    a_tags.makeEmpty();

    // Fill velocity ghosts.
    {
        if (s_nu > 0.0) {
            // Use viscous BCs
            VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
            this->setGhostsVelocity(*m_vel_new_ptr, velBC, m_time);
        } else {
            // Use inviscid BCs for tracing
            VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
            this->setGhostsVelocity(*m_vel_new_ptr, velBC, m_time);
        }
        // TODO: Use copier.
        m_vel_new_ptr->exchange();
    }


    // ------ Begin of specialized tagging ------
    // To eliminate the painful removal of //s in emacs and vi, I've
    // put all of these in if(0) / if(1) statements.

    // Random tagging - for debugging.
    if (0) {
        // Should we tag on this level?
        int levelTaggingLuck = rand() % 10 + 1;
        pout() << "levelTaggingLuck = " << levelTaggingLuck << endl;
        if (levelTaggingLuck > 3) {
            const IntVect& smallEnd = m_problem_domain.domainBox().smallEnd();
            const IntVect& size = m_problem_domain.size();

            const int maxLevelTags = 4*numProc();//(m_level == 0 ? numProc(): 1);
            for (int idx = 0; idx < maxLevelTags; ++idx) {
                IntVect levelTag;
                for (int dir = 0; dir < SpaceDim; ++dir) {
                    levelTag[dir] = rand() % size[dir] + smallEnd[dir];
                }
                a_tags |= levelTag;
            }
        }

        // if (s_tags_grow > 0) {
        //     a_tags.grow(s_tags_grow);
        // }

        return;
    }

    // L-shaped tag in middle of domain
    if (0) {
        Box region;

        for (int dir = 0; dir < SpaceDim; ++dir) {
            region = m_problem_domain.domainBox();

            for (int growDir = 0; growDir < SpaceDim; ++growDir) {
                region.growLo(growDir, -m_problem_domain.size(growDir) / 2.6);
                if (growDir != dir) {
                    region.growHi(growDir, -m_problem_domain.size(growDir) / 2.6);
                }
            }

            a_tags |= region;
        }

        return;
    }

    // Tag top and bottom boundaries
    if (0) {
        // Box region = m_problem_domain.domainBox();
        // D_TERM(,
        //        region.setSmall(0, region.bigEnd(0) - region.size(0)/2);,
        //        region.setSmall(1, region.bigEnd(1) - region.size(1)/2);)
        // region = bdryHi(region, SpaceDim-1, 1);
        // region.shiftHalf(SpaceDim-1, -1);

        // // region.grow(1, -region.size(0)/3);

        // a_tags |= region;

        int dir = SpaceDim-1;
        Box region = bdryLo(m_problem_domain, dir, 1);
        region.shiftHalf(dir, 1);
        a_tags |= region;

        return;
    }

    // Tag top and right boundaries
    if (0) {
        int dir = SpaceDim-1;
        Box region = bdryHi(m_problem_domain, dir, 1);
        region.shiftHalf(dir, -1);
        a_tags |= region;

        // if (m_level == 1) {
        //     dir = 0;
        //     region = bdryHi(m_problem_domain, dir, 1);
        //     region.shiftHalf(dir, -1);
        //     a_tags |= region;
        // }

        return;
    }

    // Tag at hill surface
    if (0) {
        // Hill surface location
        const Box domBox = m_problem_domain.domainBox();

        IntVect hillTop = domBox.smallEnd();
        hillTop += domBox.size() / 2;
        hillTop[SpaceDim-1] = domBox.smallEnd(SpaceDim-1);

        // Initialize lo and hi points at center of hill top.
        IntVect loEnd = hillTop;
        IntVect hiEnd = hillTop;

        // Streamwise extents.
        loEnd[0] = hillTop[0] - domBox.size(0) / 20;
        hiEnd[0] = hillTop[0] + domBox.size(0) / 20;

        // Vertical extents.
        hiEnd[SpaceDim-1] = hillTop[SpaceDim-1] + domBox.size(SpaceDim-1) / 4;

        // Spanwise extents.
        if (SpaceDim == 3 && m_problem_domain.isPeriodic(1)) {
            loEnd[1] = domBox.smallEnd(1);
            hiEnd[1] = domBox.bigEnd(1);
        }

        Box tagBox(loEnd, hiEnd);
        a_tags |= tagBox;

        // return;
    }
    // ------- End of specialized tagging -------


    // Tag on Richardson number
    if (s_do_Ri_tagging) {
        LevelData<FArrayBox> Ri(grids, 1);
        this->computeRiNumber(Ri, 0, m_time);

        DataIterator dit = grids.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            const FArrayBox& RiFAB = Ri[dit];
            const Box& valid = grids[dit];

            // Tag where Ri <= tol.
            BoxIterator bit(valid);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();
                if (RiFAB(cc) <= s_Ri_tag_tol) {
                    a_tags |= cc;
                }
            } // end loop over valid (bit)
        } // end loop over grids (dit)
    } // end if tagging on Gradient Richardson number


    // Tag on vorticity
    // We will tag where ever |vort| / max|vort| >= s_magvort_quota
    // Note that there will always be some portion of the domain that gets
    // tagged -- that's why this is a quota.
    if (s_magvort_tag_quota > 0.0) {
        DataIterator dit = grids.dataIterator();

#       if CH_SPACEDIM == 2
            // Compute vorticity
            LevelData<FArrayBox> magVort(grids, 1);
            this->computeVorticity(magVort);

            // Compute |vorticity|
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox& magVortFAB = magVort[dit];
                const FArrayBox& gdnFAB = m_levGeoPtr->getCCgdn()[dit];
                const FArrayBox& JinvFAB = m_levGeoPtr->getCCJinv()[dit];
                const Box region = magVortFAB.box();
                const int gdnComp = LevelGeometry::symTensorCompCC(SpaceDim-1, SpaceDim-1);

                FORT_TWOFORMMAG2D(
                    CHF_FRA1(magVortFAB,0),
                    CHF_CONST_FRA1(magVortFAB,0),
                    CHF_CONST_FRA1(gdnFAB,gdnComp),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(region));
            }

#       elif CH_SPACEDIM == 3
            // Compute vorticity
            LevelData<FArrayBox> vorticity(grids, 3);
            this->computeVorticity(vorticity);

            // Compute |vorticity|
            LevelData<FArrayBox> magVort(grids, 1);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox& magVortFAB = magVort[dit];
                const FArrayBox& vortFAB = vorticity[dit];
                const Box region = magVortFAB.box();
                const Real p = 0.5;

                m_levGeoPtr->contractVectors(magVortFAB,
                                             vortFAB,
                                             vortFAB,
                                             dit());

                FORT_POWFAB(CHF_FRA(magVortFAB),
                            CHF_BOX(region),
                            CHF_CONST_REAL(p));
            }
#       else
#           Bad SpaceDim
#       endif

        // Calculate tagLevel = s_magvort_tag_quota * max|vort|
        Real tagLevel = computeUnmappedNorm(magVort, NULL, *m_levGeoPtr, 0);
        tagLevel *= s_magvort_tag_quota; // TODO: Should we divide by dx or dA?
        tagLevel = Abs(tagLevel);

        // Test each cell for |vort| >= tagLevel
        if (tagLevel > 0.0) {
            DataIterator dit = magVort.dataIterator();
            for (dit.reset(); dit.ok(); ++dit) {
                const FArrayBox& magVortFAB = magVort[dit];
                BoxIterator bit(magVortFAB.box());

                for (bit.begin(); bit.ok(); ++bit) {
                    const IntVect& iv = bit();

                    if (abs(magVortFAB(iv)) >= tagLevel)
                        a_tags |= iv;
                }
            } // end loop over grids
        } // end if tagLevel > 0.0
    } // end if tagging on vorticity


    // New vorticity tagging strategy. Tagging occurs where the area weighted,
    // absolute vorticity is above some threshold. This version is not a quota!
    // Also note that in 2D, we only use the "z-component" -- s_vort_tag_tol[2].
    if (s_vort_tag_tol[0] + s_vort_tag_tol[1] + s_vort_tag_tol[2] > 0.0) {
        DataIterator dit = grids.dataIterator();
        const RealVect& dx = m_levGeoPtr->getDx();
        const int vortComps = D_TERM(0,+1,+2);

        // Compute vorticity
        LevelData<FArrayBox> vort(grids, vortComps);
        this->computeVorticity(vort);

        for (dit.reset(); dit.ok(); ++dit) {
            BoxIterator bit(grids[dit]);

            if (CH_SPACEDIM == 2) {
                // Scale with dA.
                vort[dit].mult(dx[0]*dx[1], 0);

                // Tag if |vort| > s_newvort_tag_tol.
                // NOTE: We use the z-component in 2D.
                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& cc = bit();

                    if (abs(vort[dit](cc,0)) >= s_vort_tag_tol[2]) {
                        a_tags |= cc;
                    }
                }

            } else {
                // Scale with dA.
                vort[dit].mult(dx[1]*dx[2], 0);
                vort[dit].mult(dx[2]*dx[0], 1);
                vort[dit].mult(dx[0]*dx[1], 2);

                // Tag if |vort| > s_newvort_tag_tol.
                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& cc = bit();

                    if (abs(vort[dit](cc,0)) >= s_vort_tag_tol[0]) {
                        a_tags |= cc;
                    } else if (abs(vort[dit](cc,1)) >= s_vort_tag_tol[1]) {
                        a_tags |= cc;
                    } else if (abs(vort[dit](cc,2)) >= s_vort_tag_tol[2]) {
                        a_tags |= cc;
                    }
                }
            }
        } // end loop over grids (dit)
    } // end if tagging with new vorticity strategy


    // Tag on velocity differences
    if (s_vel_tag_tol > 0.0) {
        const LevelData<FArrayBox>& velocity = *m_vel_new_ptr;
        const RealVect& dx = m_levGeoPtr->getDx();

        LevelData<FArrayBox>& velRef = (LevelData<FArrayBox>&)velocity;
        LevelData<FluxBox>& JgupRef = (LevelData<FluxBox>&)m_levGeoPtr->getFCJgup();
        if (s_nu > 0.0) {
            // Use viscous BCs
            VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
            velBC.setGhosts(velRef,         // state
                            NULL,           // extrap
                            dx,             // dx
                            &JgupRef,       // JgupPtr
                            false);         // inhomogeneous
        } else {
            // Use inviscid BCs for tracing
            VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
            velBC.setGhosts(velRef,         // state
                            NULL,           // extrap
                            dx,             // dx
                            &JgupRef,       // JgupPtr
                            false);         // inhomogeneous
        }
        // velRef.exchange();

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            const Box region = grids[dit];
            const FArrayBox& velFAB = velocity[dit];

            for (int idir = 0; idir < SpaceDim; ++idir) {
                // The tags are made at the faces adjoining two cells.
                // Alter region so it is essentially FC.
                Box tagRegion = region;
                tagRegion.growHi(idir, 1);

                Real D_DECL(dif0, dif1, dif2);
                bool doTag;

                for (BoxIterator bit(tagRegion); bit.ok(); ++bit) {
                    const IntVect& right = bit();
                    const IntVect  left  = bit() - BASISV(idir);

                    D_TERM(dif0 = abs(velFAB(right,0) - velFAB(left,0));,
                           dif1 = abs(velFAB(right,1) - velFAB(left,1));,
                           dif2 = abs(velFAB(right,2) - velFAB(left,2));)

                    D_TERM(doTag = (abs(dif0) >= s_vel_tag_tol),
                                || (abs(dif1) >= s_vel_tag_tol),
                                || (abs(dif2) >= s_vel_tag_tol));

                    if (doTag) {
                        a_tags |= left;
                        a_tags |= right;
                    }
                } // end loop over box
            } // end loop over difference dirs
        } // end loop over grids
    } // end if tagging on velocity differences

    // Tag on buoyancy differences
    if (s_buoyancy_tag_tol > 0.0 && s_num_scal_comps > 0) {
        CH_assert(m_scal_new[0] != NULL);

        LevelData<FArrayBox> buoyancy;
        this->fillScalars(buoyancy, m_time, 0);
        // buoyancy.exchange();

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            const Box region = grids[dit];
            const FArrayBox& bFAB = buoyancy[dit];

            for (int idir = 0; idir < SpaceDim; ++idir) {
                // The tags are made at the faces adjoining two cells.
                // Alter region so it is "essentially" FC.
                Box tagRegion = region;
                tagRegion.growHi(idir, 1);

                Real dif;

                for (BoxIterator bit(tagRegion); bit.ok(); ++bit) {
                    const IntVect& right = bit();
                    const IntVect  left  = bit() - BASISV(idir);

                    dif = abs(bFAB(right) - bFAB(left));

                    if (dif >= s_buoyancy_tag_tol) {
                        a_tags |= left;
                        a_tags |= right;
                    }
                } // end loop over box
            } // end loop over difference dirs
        } // end loop over grids
    } // end if tagging on buoyancy differences

    // This tags via the buoyancy projected onto the vertical structure.
    TODO(); // Put s_bphi0_tag_tol into ProblemContext, etc...
    if (0) {
        static Real s_bphi0_tag_tol = 0.1;
        DataIterator dit = grids.dataIterator();

        // Get the buoyancy and vertical structure function.
        LevelData<FArrayBox> buoyancy;
        this->fillScalars(buoyancy, m_time, 0);

        LevelData<FArrayBox> phi0(grids, 1);
        this->fillVerticalStructure(phi0);

        // Project the buoyancy perturbation onto phi0.
        for (dit.reset(); dit.ok(); ++dit) {
            const FArrayBox& phiFAB = phi0[dit];
            const FArrayBox& bFAB = buoyancy[dit];
            const FArrayBox& JinvFAB = m_levGeoPtr->getCCJinv()[dit];
            const Box& valid = grids[dit];
            const int hiZ = valid.bigEnd(SpaceDim-1);

            const Real intScale = 1.0 / Real(valid.size(SpaceDim-1));

            Box bottomBox = valid;
            bottomBox.setBig(SpaceDim-1, bottomBox.smallEnd(SpaceDim-1));

            FArrayBox projFAB(valid, 1);
            projFAB.copy(bFAB, valid);
            projFAB.divide(JinvFAB, valid, 0, 0, 1);
            projFAB.mult(phiFAB, valid, 0, 0, 1);

            BoxIterator bit(bottomBox);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& startIV = bit();

                Box intBox(startIV, startIV);
                intBox.setBig(SpaceDim-1, hiZ);

                Real projVal = projFAB.sum(intBox, 0, 1) * intScale;
                if (abs(projVal) > s_bphi0_tag_tol) {
                    a_tags |= intBox;
                }
                // projFAB.setVal(projVal, intBox, 0, 1);
            }
        }
    }

    // Tag on pressure differences
    if (s_pressure_tag_tol > 0.0) {
        // Change this to false as soon as one reason arises.
        bool doPresTagging = true;

        // Do nothing if incompressible
        if (doPresTagging) {
            if (!s_isIncompressible) {
                MayDay::Warning("Cannot tag on pressure for a compressible flow");
                s_pressure_tag_tol = 0.0;
                doPresTagging = false;
            }
        }

        // Has the pressure been properly defined?
        // I figure if the pressure data holder is defined and compatible with
        // this level's grids, then we have initialized the pressure field. The
        // worst case scenario is that pres = 0 everywhere and nothing gets
        // tagged. At least it's not a seg fault.
        if (doPresTagging) {
            if (!m_ccPressure.isDefined()) {
                doPresTagging = false;
            }
        }
        if (doPresTagging) {
            if (!m_ccPressure.getBoxes().compatible(grids)) {
                doPresTagging = false;
            }
        }

        if (doPresTagging) {
            // iterate on the grids
            for (DataIterator dit(grids); dit.ok(); ++dit) {
                const FArrayBox& presFAB = m_ccPressure[dit];
                const Box& region = grids[dit];

                // iterate on the directions
                for (int idir = 0; idir < SpaceDim; ++idir) {
                    Box tagRegion = region;
                    tagRegion.growHi(idir,1);

                    Real dif;
                    for (BoxIterator bit(tagRegion); bit.ok(); ++bit) {
                        const IntVect& right = bit();
                        const IntVect  left  = bit() - BASISV(idir); // note this one not a pointer, since we subtract the basis vector

                        dif = abs(presFAB(right)-presFAB(left));

                        if (dif >= s_pressure_tag_tol) {
                            a_tags |= left;
                            a_tags |= right;
                        }
                    } // end loop over box
                } // end loop over difference dirs
            } // end loop over grids
        } // end if pressure is well-defined
    } // End tagging on pressure differences


    // Grow tags if needed
    if (s_tags_grow > 0) {
        a_tags.grow(s_tags_grow);
    }

    // This fixes the problem of extrapolating ghosts over a periodic dir.
    // It would be better to fix the extrapolator.
    // Unfortunately, this is a serial operation, but it only occurs over
    // a 2D set of cells.
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (!m_problem_domain.isPeriodic(dir)) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            const int isign = sign(iside);
            const int newPos = m_problem_domain.domainBox().sideEnd(flip(iside))[dir];

            // Get the boundary box that we will pull tagged cells from.
            Box ccBdry = bdryBox(m_problem_domain.domainBox(), dir, iside, 1);
            ccBdry.shiftHalf(dir, -isign);

            // Find tagged cells and mimic on other side.
            IntVectSet oldBdryTags = ccBdry & a_tags;
            IVSIterator ivsit(oldBdryTags);
            for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                IntVect iv = ivsit();
                iv[dir] = newPos;
                a_tags |= iv;
            } // end loop over tagged cells (ivsit)
        } // end loop over sides (sit)
    } // end loop over directions (dir)


    // Vertically extrude tags
    if (s_vert_extrude_tags) {
        const int loZ = m_problem_domain.domainBox().smallEnd(SpaceDim-1);
        const int hiZ = m_problem_domain.domainBox().bigEnd  (SpaceDim-1);

        IntVectSet extrudedTags = a_tags;

        IVSIterator ivsit(extrudedTags);
        for (ivsit.reset(); ivsit.ok(); ++ivsit) {
            const IntVect& cc = ivsit();

            Box vertStrip(cc,cc);
            vertStrip.shift(SpaceDim-1, loZ-cc[SpaceDim-1]);
            vertStrip.setBig(SpaceDim-1, hiZ);

            a_tags |= vertStrip;
        }
    }

    // Extrude tags in periodic y-dir
    if (0) {
        if (SpaceDim == 3 && m_problem_domain.isPeriodic(1)) {
            const int loY = m_problem_domain.domainBox().smallEnd(1);
            const int hiY = m_problem_domain.domainBox().bigEnd  (1);

            IntVectSet extrudedTags = a_tags;

            IVSIterator ivsit(extrudedTags);
            for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                const IntVect& cc = ivsit();

                Box yStrip(cc,cc);
                yStrip.shift(1, loY-cc[1]);
                yStrip.setBig(1, hiY);

                a_tags |= yStrip;
            }
        }
    }



    // LEAVE ME ALONE.
    // An error will be thrown if a_tags internally represented by a
    // TreeIntVectSet!!!! These lines should convert its storage to a
    // DenseIntVectSet.
    TODO();
    // a_tags |= m_problem_domain.domainBox().bigEnd() + 10*IntVect::Unit;
    a_tags &= m_problem_domain;
}


// -----------------------------------------------------------------------------
// Perform any pre-regridding ops -- lBase is the finest unchanged level.
// This function will be called from the finest level downward to a_lBase.
// NOTE: If a new level is about to be created, this function will be run on
// that level as well...even if it doesn't exist yet!
// -----------------------------------------------------------------------------
void AMRNavierStokes::preRegrid (int a_lBase, const Vector<Vector<Box> >& a_newGrids)
{
    // Do nothing.
}


// -----------------------------------------------------------------------------
// Set up data on this level after regridding.
// Called from coarsest regridded level --> finest regridded level
// Load balancing occurs here.
// -----------------------------------------------------------------------------
void AMRNavierStokes::regrid (const Vector<Box>& a_new_grids)
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::regrid " << m_level
               << " Nboxes=" << a_new_grids.size() << endl;
    }
    CH_TIME("AMRNavierStokes::regrid");

    if (s_verbosity >= 4) {
        pout() << "new grids = " << a_new_grids << endl;
    }

    // Copy this level's new grids
    m_level_grids = a_new_grids;

    // From Wikipedia: Morton ordering maps multidimensional data to one
    // dimension while preserving locality of the data points.
    mortonOrdering(m_level_grids);

    // Now balance the load
    const DisjointBoxLayout grids = loadBalance(m_level_grids);

    if (m_level_grids.size() > 0) {
        // This level is not empty...
        // Set flag to indicate that grids exist on this level
        m_is_empty = false;

        // May be necessary to set finestLevel flag
        // the way this works is that if the finer level is empty, then
        // we define this level as the finest level; if it turns out that we
        // haven't gotten to that level yet, but it will have grids, then we
        // correct this in the finer level's regrid call
        if (fineNSPtr() != NULL) {
            if (fineNSPtr()->isEmpty()) {
                // This is the finest extant level
                finestLevel(true);
            } else {
                // This is NOT the finest extant level
                finestLevel(false);
            }
        } else {
            // This is the finest extant level
            finestLevel(true);
        }

        // No matter what, if we get here, the coarser level is not the finest.
        crseNSPtr()->finestLevel(false);

        // Removed Copier::define code from here. 2016-3-16 ES

        // Setup post-regrid smoothers.
        if (s_smooth_after_regrid) {
            AMRNavierStokes* baseNSPtr = this;
            if (m_level > 0) baseNSPtr = baseNSPtr->crseNSPtr();

            if (baseNSPtr->m_regrid_smoothing_done == false) {
                baseNSPtr->setupPostRegridSmoothing(baseNSPtr->m_level);
            }
        }

        // since we're using pointers here, it's easy to save old
        // data -- then clean up afterwards
        LevelData<FArrayBox>* old_newVelPtr = m_vel_new_ptr;
        LevelData<FArrayBox>* old_oldVelPtr = m_vel_old_ptr;
        LevelData<FArrayBox>* old_newLambdaPtr = m_lambda_new_ptr;
        LevelData<FArrayBox>* old_oldLambdaPtr = m_lambda_old_ptr;
        Vector<LevelData<FArrayBox>*> old_oldScals(m_scal_old.size(),NULL);
        Vector<LevelData<FArrayBox>*> old_newScals(m_scal_new.size(),NULL);

        // Loop over scalar components
        for (int comp = 0; comp < m_scal_old.size(); ++comp) {
            old_oldScals[comp] = m_scal_old[comp];
            old_newScals[comp] = m_scal_new[comp];
        }

        // reshape state with new grids
        IntVect ghostVect(D_DECL(1,1,1));
        m_vel_new_ptr = new LevelData<FArrayBox>(grids, CH_SPACEDIM, ghostVect);
        m_vel_old_ptr = new LevelData<FArrayBox>(grids, CH_SPACEDIM, ghostVect);

        m_lambda_new_ptr = new LevelData<FArrayBox>(grids, 1, ghostVect);
        m_lambda_old_ptr = new LevelData<FArrayBox>(grids, 1, ghostVect);

        m_scal_new.resize(s_num_scal_comps);
        m_scal_old.resize(s_num_scal_comps);

        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            m_scal_new[comp] = new LevelData<FArrayBox>(grids, 1, ghostVect);
            m_scal_old[comp] = new LevelData<FArrayBox>(grids, 1, ghostVect);
        }

        m_macPressure.define(grids, 1, IntVect::Unit);
        m_ccPressure.define(grids, 1, IntVect::Unit);
        m_syncPressure.define(grids, 1, IntVect::Unit);
        m_eLambda.define(grids, 1, IntVect::Unit);
        m_gradELambda.define(grids, 1, IntVect::Unit);

        m_ccPressureState = CCPressureState::BOGUS;
        m_syncPressureState = SyncPressureState::BOGUS;
        m_eLambdaState = ELambdaState::BOGUS;
        m_gradELambdaState = GradELambdaState::BOGUS;

        this->setAllBogus();

        // Set up data structures
        this->levelSetup(grids);

        // interpolate from coarser level
        if (m_coarser_level_ptr != NULL) {
            AMRNavierStokes* amrns_ptr = this->crseNSPtr();
            const IntVect& nRefCrse = amrns_ptr->m_ref_ratio;
            const bool considerCellVols = true;

            { // velocity
                MappedFineInterp fine_interp(grids, CH_SPACEDIM, nRefCrse, m_problem_domain,
                                             m_levGeoPtr,
                                             considerCellVols);
                fine_interp.interpToFine(newVel(), amrns_ptr->newVel());
                fine_interp.interpToFine(oldVel(), amrns_ptr->oldVel());
            }

            { // lambda and scalars
                MappedFineInterp fine_interp_scal(grids, 1, nRefCrse, m_problem_domain,
                                                  m_levGeoPtr,
                                                  considerCellVols);

                fine_interp_scal.interpToFine(newLambda(), amrns_ptr->newLambda());

                for (int comp = 0; comp < s_num_scal_comps; ++comp) {
                    fine_interp_scal.interpToFine(newScal(comp), amrns_ptr->newScal(comp));
                    fine_interp_scal.interpToFine(oldScal(comp), amrns_ptr->oldScal(comp));
                }
            }
        } // end if there is a coarser level

        // copy from old state
        if (old_newVelPtr != NULL) {
            old_newVelPtr->copyTo(old_newVelPtr->interval(),
                                  newVel(),
                                  newVel().interval());
        }

        if (old_oldVelPtr != NULL) {
            old_oldVelPtr->copyTo(old_oldVelPtr->interval(),
                                  oldVel(),
                                  oldVel().interval());
        }

        if (old_newLambdaPtr != NULL) {
            old_newLambdaPtr->copyTo(old_newLambdaPtr->interval(),
                                     newLambda(),
                                     newLambda().interval());
        }

        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            if (old_oldScals[comp] != NULL) {
                old_oldScals[comp]->copyTo(old_oldScals[comp]->interval(),
                                           *m_scal_old[comp],
                                           m_scal_old[comp]->interval());
            }

            if (old_newScals[comp] != NULL) {
                old_newScals[comp]->copyTo(old_newScals[comp]->interval(),
                                           *m_scal_new[comp],
                                           m_scal_new[comp]->interval());
            }
        } // end loop over scalar components

        // clean up before these pointers go out of scope...
        if (old_newVelPtr != 0) {
            delete old_newVelPtr;
            old_newVelPtr = NULL;
        }

        if (old_oldVelPtr != NULL) {
            delete old_oldVelPtr;
            old_oldVelPtr = NULL;
        }

        if (old_newLambdaPtr != NULL) {
            delete old_newLambdaPtr;
            old_newLambdaPtr = NULL;
        }

        if (old_oldLambdaPtr != NULL) {
            delete old_oldLambdaPtr;
            old_oldLambdaPtr = NULL;
        }

        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            if (old_oldScals[comp] != NULL) {
                delete old_oldScals[comp];
                old_oldScals[comp]  = NULL;
            }

            if (old_newScals[comp] != NULL) {
                delete old_newScals[comp];
                old_newScals[comp] = NULL;
            }
        } // end loop over components for scalars

    } else {
        // The new level is empty -- just clear everything

        // This level is empty...
        // Set flag to indicate that grids don't exist on this level
        m_is_empty = true;

        // should just be able to delete pointers here
        if (m_vel_new_ptr != NULL) {
            delete m_vel_new_ptr;
            m_vel_new_ptr = NULL;
        }

        if (m_vel_old_ptr != NULL) {
            delete m_vel_old_ptr;
            m_vel_old_ptr = NULL;
        }

        if (m_lambda_new_ptr != NULL) {
            delete m_lambda_new_ptr;
            m_lambda_new_ptr = NULL;
        }

        if (m_lambda_old_ptr != NULL) {
            delete m_lambda_old_ptr;
            m_lambda_old_ptr = NULL;
        }

        for (int comp = 0; comp < m_scal_new.size(); ++comp) {
            if (m_scal_new[comp] != NULL) {
                delete m_scal_new[comp];
                m_scal_new[comp] = NULL;
            }

            if (m_scal_old[comp] != NULL) {
                delete m_scal_old[comp];
                m_scal_old[comp] = NULL;
            }
        } // end loop over scalar components

        m_macPressure.clear();
        m_ccPressure.clear();
        m_syncPressure.clear();
        m_eLambda.clear();
        m_gradELambda.clear();

        m_ccPressureState = CCPressureState::UNDEFINED;
        m_syncPressureState = SyncPressureState::UNDEFINED;
        m_eLambdaState = ELambdaState::UNDEFINED;
        m_gradELambdaState = GradELambdaState::UNDEFINED;

        // This also implies that this is no longer the finest level
        this->finestLevel(false);

        // This was not in the orig code.
        // std::cout << "Severing the link between this levGeo and the coarser." << endl;
        m_levGeoPtr->setFinerPtr(NULL);
        m_levGeoPtr->setCoarserPtr(NULL);

        // if a coarser level exists, let it know that it's the finest level
        if (m_coarser_level_ptr != NULL) {
            AMRNavierStokes* crse_amrns_ptr = dynamic_cast<AMRNavierStokes*>(m_coarser_level_ptr);
            if (crse_amrns_ptr != NULL) {
                if (!crse_amrns_ptr->isEmpty()) {
                    crse_amrns_ptr->finestLevel(true);
                }
            } else {
                MayDay::Error("in AMRNavierStokes::regrid: m_coarser_level_ptr is not castable to AMRNavierStokes*");
            }
        }
    } // end if new level is empty
}


// -----------------------------------------------------------------------------
// Performs post-regrid smoothing, projection, and pressure initialization.
// a_lBase is the finest unchanged level.
// Called from finest regridded level --> a_lBase (coarsest regridded level - 1)
// -----------------------------------------------------------------------------
void AMRNavierStokes::postRegrid (int a_lBase)
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::postRegrid on level " << m_level
               << " with a_lBase = " << a_lBase << endl;
    }
    CH_TIME("AMRNavierStokes::postRegrid");

    // This looks remarkably similar to postInitialize.
    // We need to re-project velocity and recompute e_lambda;

    // Do this from the finest extant level
    if (finestLevel() && (m_level != a_lBase)) {
        // This seems to make no difference.
        // m_levGeoPtr->averageMetricsDown();

        const int numLevels = m_level + 1;

        // Find a pointer to lbase
        AMRNavierStokes* thisLevelData = this;
        for (signed int lev = m_level; lev > a_lBase; --lev) {
            thisLevelData = thisLevelData->crseNSPtr();
        }
        CH_assert(thisLevelData->m_level == a_lBase);
        AMRNavierStokes *const lBaseAMRNavierStokes = thisLevelData;


        // 1. Smoothing...

        // Smooth the composite velocity field (from lbase up) if needed
        if (s_smooth_after_regrid) {
            lBaseAMRNavierStokes->doPostRegridSmoothing(a_lBase);
        }


        // 2. Projection...

        // Collect velocity pointers from lbase, up.
        Vector<LevelData<FArrayBox>*> amrVel = gatherNewVelStartingWith(lBaseAMRNavierStokes);
        CH_assert(amrVel.size() == numLevels);

        // We may need the coarser CFBCs -- The coarse data will not
        // need ghost cells if base level is properly nested.
        LevelData<FArrayBox> crseVel;
        if (a_lBase > 0) {
            // Get coarse grids
            AMRNavierStokes* crseLevelPtr = lBaseAMRNavierStokes->crseNSPtr();
            const DisjointBoxLayout& crseGrids = crseLevelPtr->newVel().getBoxes();

            // Fill the coarse velocity (fillVelocity does not define the data holder..)
            crseVel.define(crseGrids, SpaceDim);
            crseLevelPtr->fillVelocity(crseVel, m_time);

            // Collect these fields in our composite containers
            amrVel[a_lBase-1] = &crseVel;
        }

        // Collect sync pressure pointers, then initialize everything to zero.
        // lbase needs to be zero as a homogeneous CF-BC. We do this without
        // clobbering the true sync pressure at lbase.
        Vector<LevelData<FArrayBox>*> eSync = lBaseAMRNavierStokes->gatherSyncPressure();
        LevelData<FArrayBox> lBase_eSync(lBaseAMRNavierStokes->newVel().getBoxes(), 1);
        eSync[a_lBase] = &lBase_eSync;
        setValLevels(eSync, a_lBase, m_level, 0.0);

        // Note that initial velocity projection is done holding
        // lBase velocity constant.
        if (s_isIncompressible && s_initial_projection_iters > 0) {
            pout() << "Init projection on levels " << a_lBase+1 << " to "
                   << m_level << ": " << std::flush;

            AMRCCProjector projObj;
            projObj.define(eSync, *m_physBCPtr, *m_levGeoPtr, NULL);

            for (int iter = 0; iter < s_initial_projection_iters; ++iter) {
                projObj.project(amrVel,
                                *m_levGeoPtr,
                                a_lBase+1,   // lmin
                                m_level,     // lmax
                                m_time,      // newTime
                                1.0,         // dt
                                false,       // velIsFlux,
                                false,       // zero-out pressure
                                false);      // force homogeneous
            }

            // m_syncPressure on this and all higher levels contains init
            // projection data that should not be used to influence dynamics.
            setSyncPressureStates(SyncPressureState::REGRID);
        }


        // 3. VD correction...
        // Compute freestream preservation correction on lbase as well
        // (since grad(eLambda) on lbase accounts for the grids on lbase+1
        lBaseAMRNavierStokes->computeVDCorrection(true);


        // 4. Initialize pressure...

        // Now that we have projected and corrected volume discrepancies (freestream
        // preservation), we can now calculate the initial pressures on the new grids
        // if needed.
        // NOTE: The Chombo folks reset the pressure on lbase as well....find out why!
        // I can say that if we don't do like Chombo, an assert will trip. Go ahead. Try it.
        if (s_initial_pressure_iters > 0) {
            lBaseAMRNavierStokes->initializeGlobalPressure();
        } else {
            if (s_isIncompressible) {
                MayDay::Warning("Fluid is incompressible and you don't want to initialize the pressure?! "
                                "You're probably going to blow up.");
            }
        }
    } // end if finestLevel && m_level > 0
}


// -----------------------------------------------------------------------------
// Constructs the RHS of the smooth operations to be performed after regridding.
// -----------------------------------------------------------------------------
void AMRNavierStokes::setupPostRegridSmoothing (int a_lBase)
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::setupPostRegridSmoothing with a_lBase = " << a_lBase << endl;
    }
    CH_TIME("AMRNavierStokes::setupPostRegridSmoothing");

    // Find the finest extant level. This is a bit tricky since we are in the
    // middle of regridding and finestLevel() cannot be trusted.
    AMRNavierStokes* thisNSPtr = this;
    while (thisNSPtr != NULL) {
        if (thisNSPtr->m_vel_new_ptr == NULL) {
            CH_assert(thisNSPtr->m_level > 0);
            thisNSPtr = thisNSPtr->crseNSPtr();
            break;

        } else if(thisNSPtr->m_vel_new_ptr->getBoxes().size() == 0) {
            CH_assert(thisNSPtr->m_level > 0);
            thisNSPtr = thisNSPtr->crseNSPtr();
            break;

        } else if(thisNSPtr->finestLevel()) {
            break;
        }

        thisNSPtr = thisNSPtr->fineNSPtr();
    }
    CH_assert(thisNSPtr != NULL);
    const int finestLevel = thisNSPtr->m_level;

    // The CF BCs will come from the level under a_lBase
    const int CFBCLev = (a_lBase > 0)? a_lBase-1: a_lBase;

    if (s_verbosity >= 4) {
        pout() << "Smoothing is being setup with a_lBase = " << a_lBase
               << " and finestLevel = " << finestLevel << endl;
    }

    // Let's get our sanity checks out of the way
#   ifndef NDEBUG
        // Should we be in this function?
        CH_assert(s_smooth_after_regrid);
        CH_assert(m_level == a_lBase);

        // Starting from CFBCLev, check all of the data structures to make
        // sure they agree with one another.
        thisNSPtr = this;
        if (a_lBase > 0) thisNSPtr = thisNSPtr->crseNSPtr();

        for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
            CH_assert(thisNSPtr != NULL);
            CH_assert(thisNSPtr->m_level == lev);
            CH_assert(thisNSPtr->m_regrid_smoothing_done == false);

            CH_assert(thisNSPtr->m_levGeoPtr != NULL);
            CH_assert(thisNSPtr->m_vel_new_ptr != NULL);

            CH_assert(thisNSPtr->m_levGeoPtr->getBoxes()        == thisNSPtr->newVel().getBoxes());
            CH_assert(thisNSPtr->m_levGeoPtr->getDomain()       == thisNSPtr->newVel().getBoxes().physDomain());
            CH_assert(thisNSPtr->m_levGeoPtr->getFineRefRatio() == thisNSPtr->refRatio()
                || (lev == finestLevel && thisNSPtr->m_levGeoPtr->getFineRefRatio() == IntVect::Unit));

            CH_assert(thisNSPtr->newVel().getBoxes().size() > 0);
            CH_assert(thisNSPtr->m_levGeoPtr->getBoxes().size() > 0);

            thisNSPtr = thisNSPtr->fineNSPtr();
        }
#   endif

    // Let's setup all of the structures that will be needed by both
    // the velocity smoother and the scalar smoothers first...

    // Declare AMR field holders
    Vector<LevelData<FArrayBox>* > amrS(finestLevel+1,NULL);
    Vector<LevelData<FArrayBox>* > amrLapS(finestLevel+1, NULL);
    Vector<MappedCoarseAverage*> amrAvgDown(finestLevel+1, NULL);

    // Loop over levels and allocate/define structures. Start at CFBCLev.
    thisNSPtr = this;
    if (a_lBase > 0) {
        thisNSPtr = thisNSPtr->crseNSPtr();
    }

    for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
        CH_assert(thisNSPtr->m_level == lev);

        // These will be aliased to individual vel comps later.
        amrS[lev] = new LevelData<FArrayBox>;
        amrLapS[lev] = new LevelData<FArrayBox>;

        // The MappedCoarseAverage utility exists on each level and sends data to the
        // _coarser_ level. This does not need to be defined to alter CFBCLev.
        if (lev > a_lBase) {
            amrAvgDown[lev] = new MappedCoarseAverage(thisNSPtr->m_levGeoPtr->getBoxes(),
                                                      thisNSPtr->crseNSPtr()->m_levGeoPtr->getBoxes(),
                                                      1,    // comps
                                                      thisNSPtr->m_levGeoPtr->getCrseRefRatio(),
                                                      IntVect::Zero);
        }

        // Move on to finer level
        thisNSPtr = thisNSPtr->fineNSPtr();
    }

    if (s_nu > 0.0) {
        // Begin velocity setup block

        // Define velocity smoothing op factory. Remember, defineRegridOp needs to know
        // what level to take dt from...which should be the finest unchanged level, a_lBase.
        MappedAMRPoissonOpFactory localPoissonOpFactory;
        this->defineRegridAMROp(localPoissonOpFactory, a_lBase, s_nu);

        // Define the AMRLevelOps.
        Vector<MappedAMRPoissonOp*> localPoissonOpPtrVec(finestLevel+1, NULL);
        thisNSPtr = this;
        for (int lev = a_lBase; lev <= finestLevel; ++lev) {
            CH_assert(thisNSPtr->m_level == lev);

            const ProblemDomain& levelDomain = thisNSPtr->newVel().getBoxes().physDomain();
            localPoissonOpPtrVec[lev] = dynamic_cast<MappedAMRPoissonOp*>(localPoissonOpFactory.AMRnewOp(levelDomain));

            // Move on to finer level
            thisNSPtr = thisNSPtr->fineNSPtr();
        }

        // For each velocity component, copy all levels of data into
        // temp storage. Then, loop over all levels and compute the
        // Laplacian of this comonent. Finally, copy Lap(vel) into
        // old_vel storage space.
        for (int dir = 0; dir < SpaceDim; ++dir) {
            Interval velComps(dir,dir);
            Interval tempComps(0,0);

            // Get fields from level with CF-BCs too.
            thisNSPtr = this;
            while (CFBCLev < thisNSPtr->m_level) {
                thisNSPtr = thisNSPtr->crseNSPtr();
            }
            CH_assert(thisNSPtr->m_level == CFBCLev);

            // Alias the current vel comp.
            for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
                CH_assert(thisNSPtr->m_level == lev);

                aliasLevelData(*amrS[lev], thisNSPtr->m_vel_new_ptr, velComps);
                aliasLevelData(*amrLapS[lev], thisNSPtr->m_vel_old_ptr, velComps);

                thisNSPtr = thisNSPtr->fineNSPtr();
            }

            // Now loop over levels and apply multilevel operator.
            // also do averaging down here if necessary
            for (int lev = a_lBase; lev <= finestLevel; ++lev) {
                LevelData<FArrayBox>& levelS = *amrS[lev];
                LevelData<FArrayBox>& levelLapS = *amrLapS[lev];

                if (lev == 0) {
                    // no coarser level
                    if (lev == finestLevel) {
                        // no finer level:  this is all there is
                        localPoissonOpPtrVec[lev]->applyOp(levelLapS, levelS);
                    } else {
                        // finer level exists
                        localPoissonOpPtrVec[lev]->AMROperatorNC(levelLapS,
                                                                 *amrS[lev + 1],
                                                                 levelS,
                                                                 false, // not homogeneous
                                                                 localPoissonOpPtrVec[lev+1]);
                    }
                } else if (lev < finestLevel) {
                    // three-level operator
                    localPoissonOpPtrVec[lev]->AMROperator(levelLapS,
                                                           *amrS[lev + 1],
                                                           levelS,
                                                           *amrS[lev - 1],
                                                           false,    // not homogeneous
                                                           localPoissonOpPtrVec[lev+1]);
                } else {
                    // no finer level
                    localPoissonOpPtrVec[lev]->AMROperatorNF(levelLapS,
                                                             levelS,
                                                             *amrS[lev - 1],
                                                             false); // not homogeneous
                }
            } // end loop over levels to apply op

            // Need to do average down from finest level on down.
            // To speed this up eventually, may want to do this
            // _after_ we copy into multicomponent velocity LDF.
            // Average new s down to invalid regions
            {
                const bool considerCellVols = true;

                thisNSPtr = this;
                while (thisNSPtr->m_level < finestLevel) {
                    thisNSPtr = thisNSPtr->fineNSPtr();
                }

                for (int lev = finestLevel; lev > a_lBase; --lev) {
                    CH_assert(amrLapS.size() > lev);
                    CH_assert(amrLapS[lev] != NULL);
                    CH_assert(amrLapS[lev-1] != NULL);
                    amrAvgDown[lev]->averageToCoarse(*amrLapS[lev-1],
                                                     *amrLapS[lev],
                                                     thisNSPtr->m_levGeoPtr,
                                                     considerCellVols);
                    thisNSPtr = thisNSPtr->crseNSPtr();
                }
            }
        } // end loop over velocity components

        // Time to clean up after ourselves
        for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
            delete amrS[lev];
            delete amrLapS[lev];
            delete localPoissonOpPtrVec[lev];
        }
    } // end velocity setup block

    // Begin scalar setup loop
    for (int scalComp = 0; scalComp < s_num_scal_comps; ++scalComp) {
        // Only do all this if scalar is diffused
        if (s_scal_coeffs[scalComp] > 0.0) {

            // Get fields from level with CF-BCs too.
            thisNSPtr = this;
            while (CFBCLev < thisNSPtr->m_level) {
                thisNSPtr = thisNSPtr->crseNSPtr();
            }
            CH_assert(thisNSPtr->m_level == CFBCLev);

            // Obtain pointers to the scalar fields
            for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
                CH_assert(thisNSPtr->m_level == lev);

                amrS[lev] = thisNSPtr->m_scal_new[scalComp];
                amrLapS[lev] = thisNSPtr->m_scal_old[scalComp];

                thisNSPtr = thisNSPtr->fineNSPtr();
            }

            // Define scalar smoothing op factory. Remember, defineRegridOp needs to know
            // what level to take dt from...which should be the finest unchanged level, a_lBase.
            MappedAMRPoissonOpFactory localPoissonOpFactory;
            this->defineRegridAMROp(localPoissonOpFactory, a_lBase, s_scal_coeffs[scalComp]);

            // Define the AMRLevelOps.
            Vector<MappedAMRPoissonOp*> localPoissonOpPtrVec(finestLevel+1, NULL);
            thisNSPtr = this;
            for (int lev = a_lBase; lev <= finestLevel; ++lev) {
                CH_assert(thisNSPtr->m_level == lev);

                const ProblemDomain& levelDomain = thisNSPtr->newVel().getBoxes().physDomain();
                localPoissonOpPtrVec[lev] = dynamic_cast<MappedAMRPoissonOp*>(localPoissonOpFactory.AMRnewOp(levelDomain));

                // Move on to finer level
                thisNSPtr = thisNSPtr->fineNSPtr();
            }

            // Now, compute laplacian and average down if necessary
            for (int lev = a_lBase; lev <= finestLevel; ++lev) {
                LevelData<FArrayBox>& levelS = *amrS[lev];
                LevelData<FArrayBox>& levelLapS = *amrLapS[lev];

                if (lev == 0) {
                    // no coarser level
                    if (lev == finestLevel) {
                        // no finer level:  this is all there is
                        localPoissonOpPtrVec[lev]->applyOp(levelLapS, levelS);
                    } else {
                        // finer level exists
                        localPoissonOpPtrVec[lev]->AMROperatorNC(levelLapS,
                                                                 *amrS[lev + 1],
                                                                 levelS,
                                                                 false, // not homogeneous
                                                                 localPoissonOpPtrVec[lev+1]);
                    }
                } else if (lev < finestLevel) {
                    // three-level operator
                    localPoissonOpPtrVec[lev]->AMROperator(levelLapS,
                                                           *amrS[lev + 1],
                                                           levelS,
                                                           *amrS[lev - 1],
                                                           false,    // not homogeneous
                                                           localPoissonOpPtrVec[lev+1]);
                } else {
                    // no finer level
                    localPoissonOpPtrVec[lev]->AMROperatorNF(levelLapS,
                                                             levelS,
                                                             *amrS[lev - 1],
                                                             false); // not homogeneous
                }
            } // end loop over levels

            // Delete ops for this scalar
            for (int lev = a_lBase; lev <= finestLevel; ++lev) {
                // Do not delete the scalars!
                delete localPoissonOpPtrVec[lev];
            }
        } // end if scalar is diffused
    } // end loop over scalar components

    // Time to clean up after ourselves
    for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
        delete amrAvgDown[lev];
    }

    // Set flags.
    thisNSPtr = this;
    while(thisNSPtr != NULL) {
        thisNSPtr->m_regrid_smoothing_done = true;
        thisNSPtr = thisNSPtr->fineNSPtr();
    }
}



// -----------------------------------------------------------------------------
// Smooths composite fields after regridding
// -----------------------------------------------------------------------------
void AMRNavierStokes::doPostRegridSmoothing (int a_lBase)
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::doPostRegridSmoothing with a_lBase = " << a_lBase << endl;
    }
    CH_TIME("AMRNavierStokes::doPostRegridSmoothing");

    // in this function, we take the fillpatched fields stored in
    // the old-time storage (which have been modified to contain
    // s - mu*lap(s) in the regrid() function) and perform an
    // elliptic solve which should result in a smoothed s

    // This function should be called on the lbase level
    CH_assert (m_level == a_lBase);

    // Find finest level
    AMRNavierStokes* thisNSPtr = this;
    while (!thisNSPtr->finestLevel()) {
        thisNSPtr = thisNSPtr->fineNSPtr();
    }
    CH_assert(thisNSPtr->finestLevel());
    CH_assert(thisNSPtr->m_levGeoPtr->getBoxes().size() > 0);
    int finestLevel = thisNSPtr->m_level;

    // Find the CF BC level
    const int CFBCLev = (a_lBase > 0)? a_lBase-1: a_lBase;

    if (s_verbosity > 3) {
        pout() << "Smoothing is being performed with a_lBase = " << a_lBase
               << " and finestLevel = " << finestLevel << endl;
    }

    // Let's get our sanity checks out of the way
#   ifndef NDEBUG
        // Make sure thisNSPtr is what we expect
        CH_assert(thisNSPtr->finestLevel());
        CH_assert(thisNSPtr->m_vel_new_ptr != NULL);
        CH_assert(thisNSPtr->m_vel_new_ptr->getBoxes().size() > 0);

        // Now, starting from CFBCLev, check all of the data structures to make
        // sure they agree with one another.
        thisNSPtr = this;
        if (a_lBase > 0) thisNSPtr = thisNSPtr->crseNSPtr();

        for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
            CH_assert(thisNSPtr != NULL);
            CH_assert(thisNSPtr->m_level == lev);
            if (CFBCLev < a_lBase && thisNSPtr->m_level == CFBCLev) {
                CH_assert(thisNSPtr->m_regrid_smoothing_done == false);
            } else {
                CH_assert(thisNSPtr->m_regrid_smoothing_done == true);
            }

            CH_assert(thisNSPtr->m_levGeoPtr != NULL);
            CH_assert(thisNSPtr->m_vel_new_ptr != NULL);

            CH_assert(thisNSPtr->m_levGeoPtr->getBoxes()        == thisNSPtr->newVel().getBoxes());
            CH_assert(thisNSPtr->m_levGeoPtr->getDomain()       == thisNSPtr->newVel().getBoxes().physDomain());
            CH_assert(thisNSPtr->m_levGeoPtr->getFineRefRatio() == thisNSPtr->refRatio()
                || (lev == finestLevel && thisNSPtr->m_levGeoPtr->getFineRefRatio() == IntVect::Unit));

            CH_assert(thisNSPtr->newVel().getBoxes().size() > 0);
            CH_assert(thisNSPtr->m_levGeoPtr->getBoxes().size() > 0);

            thisNSPtr = thisNSPtr->fineNSPtr();
        }

        // If we get here without tripping an assert, then all of the levels
        // from CFBCLev to finestLevel are well defined with non-empty grids.
#   endif

    // Find coarsest ProblemDomain (needed by AMRMultiGrid)
    ProblemDomain* crseDomainPtr = NULL;
    thisNSPtr = this;
    while (thisNSPtr->m_level != 0) {
        thisNSPtr = thisNSPtr->crseNSPtr();
    }
    CH_assert(thisNSPtr->m_level == 0);
    crseDomainPtr = &(thisNSPtr->m_problem_domain);

    // Declare AMR field holders
    Vector<LevelData<FArrayBox>*> oldS(finestLevel+1, NULL);
    Vector<LevelData<FArrayBox>*> newS(finestLevel+1, NULL);
    Vector<MappedCoarseAverage*> amrAvgDown(finestLevel+1,NULL);

    // Loop over levels and allocate/define structures. Start at CFBCLev.
    thisNSPtr = this;
    if (a_lBase > 0) {
        thisNSPtr = thisNSPtr->crseNSPtr();
    }

    for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
        CH_assert(thisNSPtr->m_level == lev);

        // These will be aliased to individual vel comps later.
        newS[lev] = new LevelData<FArrayBox>;
        oldS[lev] = new LevelData<FArrayBox>;

        // The MappedCoarseAverage utility exists on each level and sends data to the
        // _coarser_ level. This does not need to be defined to alter CFBCLev.
        if (lev > a_lBase) {
            amrAvgDown[lev] = new MappedCoarseAverage(thisNSPtr->m_levGeoPtr->getBoxes(),
                                                      1,    // comps
                                                      thisNSPtr->m_levGeoPtr->getCrseRefRatio());
        }

        // Move on to finer level
        thisNSPtr = thisNSPtr->fineNSPtr();
    }

    if (s_nu > 0) {
        // Begin velocity smoothing block...

        // Define velocity smoothing op factory. Remember, defineRegridOp needs to know
        // what level to take dt from...which should be the finest unchanged level, a_lBase.
        MappedAMRPoissonOpFactory localPoissonOpFactory;
        defineRegridAMROp(localPoissonOpFactory, a_lBase, s_nu);

        // Set up bottom solver for AMRMultigrid.
        BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
        bottomSolver.m_eps = s_viscous_bottom_eps;
        bottomSolver.m_reps = s_viscous_bottom_reps;
        bottomSolver.m_imax = s_viscous_bottom_imax;
        bottomSolver.m_numRestarts = s_viscous_bottom_numRestarts;
        bottomSolver.m_hang = s_viscous_bottom_hang;
        bottomSolver.m_small = s_viscous_bottom_small;
        bottomSolver.m_verbosity = s_viscous_bottom_verbosity;
        bottomSolver.m_normType = s_viscous_bottom_normType;

        // Set up the AMRMultigrid solver.
        MappedAMRMultiGrid<LevelData<FArrayBox> > streamSolver;
        streamSolver.define(*crseDomainPtr,
                            localPoissonOpFactory,
                            &bottomSolver,
                            finestLevel + 1);
        streamSolver.m_verbosity = s_viscous_AMRMG_verbosity;
        streamSolver.m_imin = s_viscous_AMRMG_imin;
        streamSolver.setSolverParameters(s_viscous_AMRMG_num_smooth_down,
                                         s_viscous_AMRMG_num_smooth_up,
                                         s_viscous_AMRMG_num_smooth_bottom,
                                         s_viscous_AMRMG_numMG,
                                         s_viscous_AMRMG_imax,
                                         s_viscous_AMRMG_eps,
                                         s_viscous_AMRMG_hang,
                                         s_viscous_AMRMG_normThresh);

        // now loop over velocity components
        for (int dir = 0; dir < SpaceDim; ++dir) {
            Interval velComps(dir,dir);
            Interval tempComps(0,0);

            // Get fields from level with CF-BCs too.
            thisNSPtr = this;
            if (CFBCLev < m_level) {
                thisNSPtr = thisNSPtr->crseNSPtr();
            }
            CH_assert(thisNSPtr->m_level == CFBCLev);

            // Alias the current vel comp.
            for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
                CH_assert(thisNSPtr->m_level == lev);
                // CH_assert(lev == CFBCLev || thisNSPtr->m_regrid_smoothing_done);

                aliasLevelData(*oldS[lev], thisNSPtr->m_vel_old_ptr, velComps);
                aliasLevelData(*newS[lev], thisNSPtr->m_vel_new_ptr, velComps);

                thisNSPtr = thisNSPtr->fineNSPtr();
            }

            // Perform exchanges.
            for (int lev = a_lBase+1; lev <= finestLevel; ++lev) {
                Copier excp;
                excp.exchangeDefine(newS[lev]->getBoxes(), newS[lev]->ghostVect());
                newS[lev]->exchange(excp);
            }

            // Create notice of the solve.
            if (s_verbosity >= 1) {
                pout() << "Smoothing velocity field from levels "
                       << std::fixed << a_lBase << " to " << finestLevel << ":"
                       << flush;
            }

            // Solve!
            streamSolver.solve(newS, oldS, finestLevel, a_lBase,
                               false); // don't initialize newS to zero

            // Average new s down to invalid regions
            {
                const bool considerCellVols = true;

                thisNSPtr = this;
                while (thisNSPtr->m_level < finestLevel) {
                    thisNSPtr = thisNSPtr->fineNSPtr();
                }

                for (int lev = finestLevel; lev > a_lBase; --lev) {
                    CH_assert(newS.size() > lev);
                    CH_assert(newS[lev] != NULL);
                    CH_assert(newS[lev-1] != NULL);

                    amrAvgDown[lev]->averageToCoarse(*newS[lev-1],
                                                     *newS[lev],
                                                     thisNSPtr->m_levGeoPtr,
                                                     considerCellVols);
                    thisNSPtr = thisNSPtr->crseNSPtr();
                }
            }
        } // end loop over velocity components

        // Since scalars won't need temp storge, clean it up here
        for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
            delete newS[lev];
            delete oldS[lev];
        }
    } // end velocity smoothing

    // now do scalars
    for (int scalComp = 0; scalComp < s_num_scal_comps; ++scalComp) {
        // only do this if scalar is diffused
        if (s_scal_coeffs[scalComp] > 0.0) {

            // Get fields from level with CF-BCs too.
            thisNSPtr = this;
            while (thisNSPtr->m_level != CFBCLev) {
                thisNSPtr = thisNSPtr->crseNSPtr();
            }
            CH_assert(thisNSPtr->m_level == CFBCLev);

            // Obtain pointers to the scalar fields
            for (int lev = CFBCLev; lev <= finestLevel; ++lev) {
                newS[lev] = thisNSPtr->m_scal_new[scalComp];
                oldS[lev] = thisNSPtr->m_scal_old[scalComp];

                thisNSPtr = thisNSPtr->fineNSPtr();
            }

            // Perform exchanges.
            for (int lev = a_lBase+1; lev <= finestLevel; ++lev) {
                Copier excp;
                excp.exchangeDefine(newS[lev]->getBoxes(), newS[lev]->ghostVect());
                newS[lev]->exchange(excp);
            }

            // Define scalar smoothing op factory. Remember, defineRegridOp needs to know
            // what level to take dt from...which should be the finest unchanged level, a_lBase.
            MappedAMRPoissonOpFactory localPoissonOpFactory;
            defineRegridAMROp(localPoissonOpFactory, a_lBase, s_scal_coeffs[scalComp]);

            // Set up bottom solver for AMRMultigrid.
            BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
            bottomSolver.m_eps = s_viscous_bottom_eps;
            bottomSolver.m_reps = s_viscous_bottom_reps;
            bottomSolver.m_imax = s_viscous_bottom_imax;
            bottomSolver.m_numRestarts = s_viscous_bottom_numRestarts;
            bottomSolver.m_hang = s_viscous_bottom_hang;
            bottomSolver.m_small = s_viscous_bottom_small;
            bottomSolver.m_verbosity = s_viscous_bottom_verbosity;
            bottomSolver.m_normType = s_viscous_bottom_normType;

            // Set up AMRMultigrid solver.
            MappedAMRMultiGrid<LevelData<FArrayBox> > streamSolver;
            streamSolver.define(*crseDomainPtr,
                                localPoissonOpFactory,
                                &bottomSolver,
                                finestLevel + 1);
            streamSolver.m_verbosity = s_viscous_AMRMG_verbosity;
            streamSolver.m_imin = s_viscous_AMRMG_imin;
            streamSolver.setSolverParameters(s_viscous_AMRMG_num_smooth_down,
                                             s_viscous_AMRMG_num_smooth_up,
                                             s_viscous_AMRMG_num_smooth_bottom,
                                             s_viscous_AMRMG_numMG,
                                             s_viscous_AMRMG_imax,
                                             s_viscous_AMRMG_eps,
                                             s_viscous_AMRMG_hang,
                                             s_viscous_AMRMG_normThresh);

            // Post notice of the solve.
            if (s_verbosity >= 1) {
                pout() << "Smoothing scalar[" << scalComp << "] field from levels "
                       << std::fixed << a_lBase << " to " << finestLevel << ":"
                       << flush;
            }

            // Solve! (Last two args are max level and base level)
            streamSolver.solve(newS, oldS, finestLevel, a_lBase,
                               false); // dont initialize newS to zero

            // Average new s down to invalid regions
            {
                const bool considerCellVols = true;

                thisNSPtr = this;
                while (thisNSPtr->m_level < finestLevel) {
                    thisNSPtr = thisNSPtr->fineNSPtr();
                }

                for (int lev = finestLevel; lev > a_lBase; --lev) {
                    CH_assert(newS.size() > lev);
                    CH_assert(newS[lev] != NULL);
                    CH_assert(newS[lev-1] != NULL);

                    amrAvgDown[lev]->averageToCoarse(*newS[lev-1],
                                                     *newS[lev],
                                                     thisNSPtr->m_levGeoPtr,
                                                     considerCellVols);
                    thisNSPtr = thisNSPtr->crseNSPtr();
                }
            }

        }  // end if this scalar is diffused
    }  // end loop over scalars

    // Time to clean up after ourselves and reset flags.
    thisNSPtr = this;
    for (int lev = a_lBase; lev <= finestLevel; ++lev) {
        delete amrAvgDown[lev];

        CH_assert(thisNSPtr->m_level == lev);
        thisNSPtr->m_regrid_smoothing_done = false;

        thisNSPtr = thisNSPtr->fineNSPtr();
    }
}


// -----------------------------------------------------------------------------
// Defines operator factory used in smoothing during regridding
// -----------------------------------------------------------------------------
void AMRNavierStokes::defineRegridAMROp (MappedAMRPoissonOpFactory&       a_factory,
                                         const int&                       a_lBase,
                                         const Real                       a_viscCoeff)
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::defineRegridAMROp on level "
               << m_level << " with a_lBase = " << a_lBase << endl;
    }

    // Want to use dt from lbase
    Real dtLBase = m_dt;
    if (a_lBase > m_level) {
        AMRNavierStokes* thisNSPtr = fineNSPtr();
        while (a_lBase > thisNSPtr->m_level) {
            thisNSPtr = thisNSPtr->fineNSPtr();
        }
        dtLBase = thisNSPtr->m_dt;
    } else if (a_lBase < m_level) {
        AMRNavierStokes* thisNSPtr = crseNSPtr();
        while (a_lBase < thisNSPtr->m_level) {
            thisNSPtr = thisNSPtr->crseNSPtr();
        }
        dtLBase = thisNSPtr->m_dt;
    }

    // define coefficient
    Real alpha = 1.0;
    Real mu = -s_regrid_smoothing_coeff * dtLBase * a_viscCoeff;
    BCMethodHolder bcHolder(m_physBCPtr->smoothingSolverBC());

    CH_assert(m_levGeoPtr != NULL);
    a_factory.define(m_levGeoPtr,
                     alpha,
                     mu,
                     bcHolder,
                     s_viscous_AMRMG_maxDepth,
                     s_viscous_AMRMG_num_smooth_precond,
                     s_viscous_AMRMG_precondMode,
                     s_viscous_AMRMG_relaxMode);
}
