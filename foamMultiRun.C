/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    foamMultiRun

Description
    Loads and executes an OpenFOAM solver modules for each region of a
    multiregion simulation e.g. for conjugate heat transfer.

    The region solvers are specified in the \c regionSolvers dictionary entry in
    \c controlDict, containing a list of pairs of region and solver names,
    e.g. for a two region case with one fluid region named
    liquid and one solid region named tubeWall:
    \verbatim
        regionSolvers
        {
            liquid          fluid;
            tubeWall        solid;
        }
    \endverbatim

    The \c regionSolvers entry is a dictionary to support name substitutions to
    simplify the specification of a single solver type for a set of
    regions, e.g.
    \verbatim
        fluidSolver     fluid;
        solidSolver     solid;

        regionSolvers
        {
            tube1             $fluidSolver;
            tubeWall1         solid;
            tube2             $fluidSolver;
            tubeWall2         solid;
            tube3             $fluidSolver;
            tubeWall3         solid;
        }
    \endverbatim

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional dual time-step coupling
    The optional \c regionTime dictionary in \c controlDict enables a
    quasi-steady gas/solid coupling.  A region with \c mode \c quasiSteady
    (typically the gas / multiphaseEuler region) is advanced with its own
    time-step in an inner loop until the interface heat flux has converged,
    while its time does NOT accumulate into the global time.  The \c transient
    regions (typically the solid) advance the global time; the number of global
    steps is therefore \c endTime / \c deltaT.

    Example:
    \verbatim
        regionTime
        {
            gas
            {
                mode            quasiSteady;
                deltaT          1e-5;
                maxIter         200;
                heatFluxTol     1e-3;
                heatFlux
                {
                    region      solid;
                    patch       solid_to_gas;
                }
            }
            solid
            {
                mode            transient;
                deltaT          5e-4;
            }
        }
    \endverbatim

Usage
    \b foamMultiRun [OPTION]

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "regionSolvers.H"
#include "pimpleMultiRegionControl.H"
#include "setDeltaT.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"
#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Per-region time control parameters read from the regionTime dictionary
struct regionTimeParams
{
    word    mode;           // "transient" (default) or "quasiSteady"
    scalar  deltaT;         // region time-step (gas pseudo-time / solid step)
    scalar  interval;       // physical time advanced per coupling step (solid)
    scalar  heatFluxTol;    // relative interface heat-flux tolerance
    label   maxIter;        // max inner iterations for the quasi-steady region
    word    heatFluxRegion; // region on which the heat flux is measured
    word    heatFluxPatch;  // patch on that region ("" = auto-detect mappedWall)
    word    heatFluxField;  // field used for the flux proxy (default "T")
};


//- Total |heat flux| proxy through a patch: sum(|snGrad(T)| * |Sf|).
//  The constant conductivity factor cancels in the relative convergence check.
scalar interfaceHeatFlux
(
    const fvMesh& mesh,
    const label patchi,
    const word& fieldName
)
{
    const volScalarField& T = mesh.lookupObject<volScalarField>(fieldName);

    scalarField qf =
        mag(T.boundaryField()[patchi].snGrad())
       *mesh.boundary()[patchi].magSf();

    return gSum(qf);
}


//- Find the region index by name (returns -1 if not found)
label regionIndex(regionSolvers& solvers, const word& name)
{
    forAll(solvers, i)
    {
        if (solvers[i].mesh.name() == name)
        {
            return i;
        }
    }

    return -1;
}


//- Perform one PIMPLE corrector iteration for the given regions using the
//  v14 solver interface (split transport predictors/correctors).
void pimpleIteration
(
    regionSolvers& solvers,
    pimpleMultiRegionControl& pimple,
    const labelList& regions
)
{
    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.flow())
        {
            s.moveMesh();
        }
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.flow())
        {
            s.motionCorrector();
        }
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.models())
        {
            s.fvModels().correct();
        }
    }

    forAll(regions, i)
    {
        solvers[regions[i]].prePredictor();
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.predictTransport() && s.pimple.flow())
        {
            s.momentumTransportPredictor();
        }
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.predictTransport() && s.pimple.thermophysics())
        {
            s.thermophysicalTransportPredictor();
        }
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.flow())
        {
            s.momentumPredictor();
        }
    }

    while (pimple.correctEnergy())
    {
        forAll(regions, i)
        {
            solver& s = solvers[regions[i]];
            if (s.pimple.thermophysics())
            {
                s.thermophysicalPredictor();
            }
        }
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.flow())
        {
            s.pressureCorrector();
        }
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.correctTransport() && s.pimple.flow())
        {
            s.momentumTransportCorrector();
        }
    }

    forAll(regions, i)
    {
        solver& s = solvers[regions[i]];
        if (s.pimple.correctTransport() && s.pimple.thermophysics())
        {
            s.thermophysicalTransportCorrector();
        }
    }
}


//- Run one quasi-steady (gas) region to heat-flux convergence without
//  accumulating its pseudo-time into the global time.
void solveGasToConvergence
(
    Time& runTime,
    regionSolvers& solvers,
    pimpleMultiRegionControl& pimple,
    const label gasRegion,
    const regionTimeParams& p
)
{
    solver& gas = solvers[gasRegion];

    // Locate the heat-flux measurement region and patch.
    // Default: measure on the gas region itself (it is the one being actively
    // solved, so its interface gradient is non-zero and converges).
    label hfRegion = regionIndex(solvers, p.heatFluxRegion);
    if (hfRegion < 0)
    {
        hfRegion = gasRegion;
    }

    if (hfRegion < 0)
    {
        FatalErrorInFunction
            << "Cannot find a region to measure the interface heat flux on"
            << exit(FatalError);
    }

    const fvMesh& hfMesh = solvers[hfRegion].mesh;

    label hfPatch = -1;
    if (p.heatFluxPatch != word::null)
    {
        hfPatch = hfMesh.boundary().findIndex(p.heatFluxPatch);
    }
    else
    {
        // Auto-detect the mapped wall patch (the CHT interface)
        forAll(hfMesh.boundary(), patchi)
        {
            if (hfMesh.boundary()[patchi].type() == "mappedWall")
            {
                hfPatch = patchi;
                break;
            }
        }
    }

    if (hfPatch < 0)
    {
        FatalErrorInFunction
            << "Cannot find the interface patch "
            << (p.heatFluxPatch == word::null ? "(auto)" : p.heatFluxPatch)
            << " on region " << hfMesh.name()
            << exit(FatalError);
    }

    // Remember the physical time only. The gas pseudo-time advances runTime
    // during the inner loop and is rewound afterwards via setTime, which
    // leaves timeIndex monotonic (this is what keeps ddt/oldTime consistent).
    const scalar tPhys = runTime.value();

    solvers.setPrefix(gasRegion);

    Info<< nl
        << "##################################################" << nl
        << "### [GAS] quasi-steady inner loop, region " << gas.mesh.name() << nl
        << "###       deltaT = " << p.deltaT << ", maxIter = " << p.maxIter
        << ", heatFluxTol = " << p.heatFluxTol << nl
        << "###       heatFlux = " << hfMesh.name() << "/"
        << hfMesh.boundary()[hfPatch].name() << "/" << p.heatFluxField << nl
        << "##################################################" << nl << endl;

    scalar qPrev = VGREAT;
    scalar qCur = 0;
    label iter = 0;

    for (iter = 0; iter < p.maxIter; ++iter)
    {
        // Advance the gas with its own pseudo-time step (global time is
        // restored at the end of the inner loop)
        runTime.setDeltaTNoAdjust(p.deltaT);

        gas.preSolve();

        runTime++;

        Info<< "    inner iteration " << iter
            << ", gas pseudo-time = " << runTime.value() << nl << endl;

        // PIMPLE corrector loop for the gas region only
        const labelList gasRegions(1, gasRegion);

        while (pimple.loop())
        {
            pimpleIteration(solvers, pimple, gasRegions);
        }

        gas.postSolve();

        // Interface heat flux after this gas step
        qCur = interfaceHeatFlux(hfMesh, hfPatch, p.heatFluxField);

        Info<< "    interface heat flux proxy = " << qCur;

        if (iter > 0)
        {
            const scalar rel =
                mag(qCur - qPrev)/(mag(qCur) + SMALL);

            Info<< ", relative change = " << rel;

            if (rel < p.heatFluxTol)
            {
                Info<< " -> converged" << nl << endl;
                ++iter;
                break;
            }
        }

        Info<< nl << endl;

        qPrev = qCur;
    }

    if (iter >= p.maxIter)
    {
        WarningInFunction
            << "Region " << gas.mesh.name()
            << " did not converge within " << p.maxIter
            << " inner iterations" << nl << endl;
    }

    // Rewind the physical time value only (timeIndex stays monotonic).
    // The gas pseudo-time is discarded; global time is NOT advanced.
    runTime.setTime(tPhys, runTime.timeIndex());

    Info<< "### [GAS] region " << gas.mesh.name()
        << " completed in " << iter << " inner iterations" << nl
        << "###       (gas pseudo-time discarded, global time NOT advanced)" << nl
        << endl;

    solvers.setGlobalPrefix();
}


//- Read the regionTime dictionary and return the per-region parameters
List<regionTimeParams> readRegionTimeParams
(
    const Time& runTime,
    regionSolvers& solvers
)
{
    List<regionTimeParams> params(solvers.size());

    const dictionary& regionTimeDict =
        runTime.controlDict().subOrEmptyDict("regionTime");

    forAll(solvers, i)
    {
        const word rname = solvers[i].mesh.name();

        params[i].mode = "transient";
        params[i].deltaT = runTime.deltaTValue();
        params[i].interval = runTime.deltaTValue();
        params[i].heatFluxTol = 1e-3;
        params[i].maxIter = 200;
        params[i].heatFluxRegion = word::null;
        params[i].heatFluxPatch = word::null;
        params[i].heatFluxField = "T";

        if (regionTimeDict.found(rname))
        {
            const dictionary& d = regionTimeDict.subDict(rname);

            params[i].mode = d.lookupOrDefault<word>("mode", "transient");
            params[i].deltaT =
                d.lookupOrDefault<scalar>("deltaT", runTime.deltaTValue());
            params[i].interval =
                d.lookupOrDefault<scalar>("interval", params[i].deltaT);
            params[i].heatFluxTol =
                d.lookupOrDefault<scalar>("heatFluxTol", 1e-3);
            params[i].maxIter =
                d.lookupOrDefault<label>("maxIter", 200);

            if (d.found("heatFlux"))
            {
                const dictionary& hf = d.subDict("heatFlux");

                params[i].heatFluxRegion =
                    hf.lookupOrDefault<word>("region", word::null);
                params[i].heatFluxPatch =
                    hf.lookupOrDefault<word>("patch", word::null);
                params[i].heatFluxField =
                    hf.lookupOrDefault<word>("field", "T");
            }
        }
    }

    return params;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // Create the region meshes and solvers
    regionSolvers solvers(runTime);

    // Create the outer PIMPLE loop and control structure
    pimpleMultiRegionControl pimple(runTime, solvers);

    // Read the optional dual time-step controls
    const List<regionTimeParams> rtParams = readRegionTimeParams(runTime, solvers);

    // Classify the regions
    labelList gasRegions;    // quasiSteady - solved to heat-flux convergence
    labelList solidRegions;  // transient  - advance the global time
    labelList allRegions;    // all region indices

    forAll(solvers, i)
    {
        allRegions.append(i);

        if (rtParams[i].mode == "quasiSteady")
        {
            gasRegions.append(i);
        }
        else
        {
            solidRegions.append(i);
        }
    }

    const bool dualMode = (gasRegions.size() > 0);

    // Number of solid sub-steps per coupling step.  For each transient (solid)
    // region the user gives a sub-step (deltaT) and a physical interval; the
    // number of sub-steps is the interval divided by deltaT, rounded UP.
    label nSolidSteps = 1;

    // Physical (solid) time-step used for the global clock.
    scalar solidDeltaT = runTime.deltaTValue();

    // Global time-step: the step of the transient (solid) regions.
    // In dual mode the gas pseudo-time never contributes to runTime.
    if (dualMode && solidRegions.size())
    {
        // Global clock = the transient (solid) region time-step.  Note: this
        // is NOT clamped by controlDict deltaT, which in dual mode only acts
        // as a default for the gas pseudo-time.
        scalar globalDeltaT = VGREAT;

        forAll(solidRegions, i)
        {
            globalDeltaT = min(globalDeltaT, rtParams[solidRegions[i]].deltaT);
        }

        solidDeltaT = globalDeltaT;
        runTime.setDeltaTNoAdjust(globalDeltaT);

        // Sub-step count from the first solid region (rounded up)
        const scalar dt = rtParams[solidRegions[0]].deltaT;
        const scalar interval = rtParams[solidRegions[0]].interval;

        nSolidSteps = max(1, label(ceil(interval/dt - 1e-6)));

        Info<< "Partitioned coupling: " << nSolidSteps
            << " solid sub-step(s) per coupling step (interval = "
            << interval << ", deltaT = " << dt << ")" << nl << endl;
    }
    else
    {
        // Standard single-clock behaviour
        setDeltaT(runTime, solvers);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    solvers.setGlobalPrefix();

    Info<< nl << "Starting time loop\n" << endl;

    label outerIter = 0;

    while (pimple.run(runTime))
    {
        Info<< nl
            << "========================================================" << nl
            << "=== OUTER COUPLING STEP " << outerIter << nl
            << "=== [SOLID] will advance global time by "
            << (nSolidSteps*runTime.deltaTValue()) << nl
            << "=== [SOLID] " << nSolidSteps << " sub-step(s) of "
            << runTime.deltaTValue() << nl
            << "========================================================" << nl << endl;

        // The list of regions solved together with the global clock.
        // In dual mode this is the solid only (the gas has already been
        // relaxed to a quasi-steady state compatible with the current solid).
        const labelList& activeRegions = dualMode ? solidRegions : allRegions;

        // --- Quasi-steady (gas) inner loops: no global-time accumulation ---
        forAll(gasRegions, gi)
        {
            const label i = gasRegions[gi];

            solveGasToConvergence(runTime, solvers, pimple, i, rtParams[i]);
        }

        // Restore the physical (solid) time-step after the gas phase, which
        // switched runTime.deltaT to the gas pseudo-time step.
        if (dualMode && solidRegions.size())
        {
            runTime.setDeltaTNoAdjust(solidDeltaT);
        }

        // --- Transient (solid) regions: advance the physical time by
        //     nSolidSteps sub-steps ---
        for (label sub = 0; sub < nSolidSteps; ++sub)
        {
            forAll(activeRegions, i)
            {
                solvers[activeRegions[i]].preSolve();
            }

            solvers.setGlobalPrefix();

            runTime++;

            Info<< nl
                << "### [SOLID] sub-step " << (sub + 1) << '/' << nSolidSteps
                << ", Time = " << runTime.userTimeName() << nl
                << "### [SOLID] solving region(s): ";

            forAll(activeRegions, i)
            {
                Info<< solvers[activeRegions[i]].mesh.name() << ' ';
            }

            Info<< nl << endl;

            // Multi-region PIMPLE corrector loop over the active regions
            while (pimple.loop())
            {
                pimpleIteration(solvers, pimple, activeRegions);
            }

            forAll(activeRegions, i)
            {
                solvers[activeRegions[i]].postSolve();
            }

            solvers.setGlobalPrefix();

            runTime.write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        ++outerIter;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
