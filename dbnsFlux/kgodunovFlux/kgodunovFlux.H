/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is not a part of foam-extend.
    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
Class
    GodunovFlux
Description
    Implementation of inviscid numerical fluxes based on Godunov scheme.
    This is a first order scheme

    Glazunov A.A., Kagenov A.M., Kostyushin K.V., Eremin I.V., Kotonogov V.A., AligasanovaK.L. (2020) 
    MATHEMATICAL MODELING OF THE INTERACTION OF A SINGLE
    SUPERSONIC JET WITH OBSTACLES. Vestnik Tomskogo gosudarstvennogo universiteta.
    Matematika i mekhanika [Tomsk State University Journal of Mathematics and Mechanics]. 63.
    pp. 87-101

    Godunov S.K., Zabrodin A.V., Ivanov M.Ya., Krayko A.N., Prokopov G.P. (1976) 
    Chislennoe reshenie mnogomernykh zadach gazovoy dinamiki [Numerical solution to multidimensional problems of gas dynamics]. Moscow: Nauka.

Author's
    K.V.Kostyushin, A.M. Kagenov
SourceFiles
    kgodunovFlux.H

\*---------------------------------------------------------------------------*/

#ifndef kgodunovFlux_H
#define kgodunovFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class godunovFlux Declaration
\*---------------------------------------------------------------------------*/
class kgodunovFlux
{

public:

    // Constructors

        //- Construct null
        kgodunovFlux()
        {}


    // Destructor - default


    // Member Functions
    //- Update fluxes based on current left and right states
    void evaluateFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& pLeft,
        const scalar& pRight,
        const vector& ULeft,
        const vector& URight,
        const scalar& TLeft,
        const scalar& TRight,
        const scalar& RLeft,
        const scalar& RRight,
        const scalar& CvLeft,
        const scalar& CvRight,
        const vector& Sf,
        const scalar& magSf
    ) const;

    void f
    (
        scalar& P, 
        scalar& p, 
        scalar& ro, 
        scalar& k, 
        scalar& f, 
        scalar& f1
    ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //