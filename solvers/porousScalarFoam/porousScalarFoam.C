/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    passiveScalarFoam

Description
    A copy of icoFoam, a transient solver for incompressible, laminar flow of Newtonian fluids. 
    This includes scalar transport of thrombin

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

       // solve(e*UEqn == -e*fvc::grad(p) - m*(sqr(e)*((1/k)*nu*U)));
	solve(s*UEqn == -s*fvc::grad(p) - m*(sqr(s)*((1/k)*nu*U)));

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {

 #	    include "ThEqn.H"
 #	    include "sEqn.H"
 #	    include "muEqn.H"
 #	    include "wallGradU.H"
 #	    include "StrainRateMag.H"

            volScalarField rAU(1.0/UEqn.A());

            U = rAU*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf())
                  + fvc::ddtPhiCorr(rAU, U, phi);


            adjustPhi(phi, U, p);
	    

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                  fvm::laplacian(rAU, p) == fvc::div(phi)	
			
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

            #include "continuityErrs.H"


            U -= rAU*fvc::grad(p);
	    
            U.correctBoundaryConditions();
        }
	
	// Solving the transport equation for thrombin
	
	solve
            (
                fvm::ddt(Th)
              + fvm::div(phi, Th)
              - fvm::laplacian(dTh, Th)
            );

	// Calculate strain rate magnitude on walls

	//forAll(wallGradU.boundaryField(), patchi)
        //     {
        //       wallGradU.boundaryField()[patchi] =
        //            -U.boundaryField()[patchi].snGrad();
        //     } 
	//	
	//     StrainRateMag = mag(wallGradU);

	// Setting a minimum threshold for thrombin of ThMin in the entire domain (avoid machine error)

	Th.max(ThMin);

	// Determine the permeability of the different cells based on thrombin concentration creating matrix of 0 and k
	// Porosity in clot region is 0.75. 0.875 halfway between 0.75 (clot region) and 1.0 (fluid region)

	s = ThCrit;	
	s = Th-s;
	s = sign(s);

	m = 0.5*(s+1.0);

	s = -0.125*s;
	s = s + 0.875;

    // Assumes entire zone is porous

	// s=0.75;
	// m=1.0;
	
        
	runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
