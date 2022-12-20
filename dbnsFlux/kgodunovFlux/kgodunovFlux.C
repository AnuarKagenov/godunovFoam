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
\*---------------------------------------------------------------------------*/

#include "kgodunovFlux.H"

//#include <stdio.h>
//#include <time.h>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kgodunovFlux::evaluateFlux
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
) const
{
    const vector normalVector = Sf / magSf;
    vector Unorm1, Unorm2;
    vector Utau1, Utau2;

    scalar W = 0.0;

    Unorm1 = (ULeft & normalVector) * normalVector;
    Unorm2 = (URight & normalVector) * normalVector;
    Utau1 = ULeft - Unorm1;
    Utau2 = URight - Unorm2;

    vector V_l = Utau1;
    vector V_p = Utau2;	
    vector Vr; 
    vector Utotal;

    scalar k = (CvRight + RRight) / CvRight;
    scalar p1 = pLeft;
    scalar p2 = pRight;
    scalar ro1 = pLeft / (RLeft * TLeft);
    scalar ro2 = pRight / (RRight * TRight);

    scalar u1 = ULeft & normalVector;
    scalar u2 = URight & normalVector;

    scalar u1_u2 = u1 - u2;
    scalar c1 = Foam::sqrt(k * p1 / ro1);
    scalar c2 = Foam::sqrt(k * p2 / ro2);

    scalar km1_2 = (k - 1.0) / 2.0;

    scalar _a1 = ro1 * c1;
    scalar _a2 = ro2 * c2;
    scalar P0 = (p1 * _a2 + p2 * _a1 + (u1_u2) * _a1 * _a2) / (_a1 + _a2);
    if ((P0 < 0.0) || (P0 > 100000000000000))
    {
        P0 = -u1_u2 + 1.0;
    }

    scalar P = P0;
    scalar f1, f2, f1_1, f2_1;

    int iter = 0;
    do
    {
        P0 = P;

        Foam::kgodunovFlux::f(P0, p1, ro1, k, f1, f1_1);
        Foam::kgodunovFlux::f(P0, p2, ro2, k, f2, f2_1);

        P = P0 - (f1 + f2 - u1_u2) / (f1_1 + f2_1);

        if(P<=0)
        {
	    P=(p1+p2)/2.0;
	    f2=0;
	    f1=0;
	    break;
        }
        
        iter++;
    } while (fabs(P - P0) > 0.000001);

    scalar U = (u1 + u2 + f2 - f1) / 2.0;

    scalar Roe = ro1; 
    scalar Ue = 0;
    scalar Pe = P; 

    if ((P >= p1) & (U > W))
    {
        scalar m1 = Foam::sqrt(ro1 / 2.0 * (P * (k + 1.0) + p1 * (k - 1.0)));
        scalar W1 = u1 - m1 / ro1;
        Vr = V_l;

        if (W > W1)
        {
            Roe = ro1 * m1 / (m1 - ro1 * (u1 - U));
            Ue = U;
            Pe = P;
        }
        else
        {
            Roe = ro1;
            Ue = u1;
            Pe = p1;
        }
    }

    if ((P >= p2) & (U < W))
    {
        Vr = V_p;

        scalar m2 = Foam::sqrt(ro2 / 2.0 * (P * (k + 1.0) + p2 * (k - 1.0)));
        scalar W2 = u2 + m2 / ro2;

        if (W < W2)
            {
                Roe = ro2 * m2 / (m2 + ro2 * (u2 - U));
                Ue = U;
                Pe = P;
            }
            else
            {
                Roe = ro2;
                Ue = u2;
                Pe = p2;
            }

    }

    if ((P < p1) & (U > W))
    {
        scalar W1 = u1 - c1;
        scalar c1_ = c1 + km1_2 * (u1 - U);
        scalar W1_ = U - c1_;
        Vr = V_l;

        if (W > W1_)
        {
            Roe = k * P / (c1_ * c1_);
            Ue = U;
            Pe = P;
        }
        else if (W < W1)
        {
            Roe = ro1;
            Ue = u1;
            Pe = p1;
        }
        else
        {
            scalar c_ = c1 * 2.0 / (k + 1.0) + (u1 - W) * (k - 1.0) / (k + 1.0);
            Roe = ro1 * Foam::pow(c_ / c1, 2.0 / (k - 1.0));
            Ue = W + c_;
            Pe = p1 * Foam::pow(c_ / c1, 2.0 * k / (k - 1.0));
        }
    }

    if ((P < p2) & (U < W))
    {
        scalar W2 = u2 + c2;
        scalar c2_ = c2 - km1_2 * (u2 - U);
        scalar W2_ = U + c2_;
        Vr = V_p;

        if (W < W2_)
        {
            Roe = k * P / (c2_ * c2_);
            Ue = U;
            Pe = P;
        }
        else if (W > W2)
        {
            Roe = ro2;
            Ue = u2;
            Pe = p2;
        }
        else
        {
            scalar c_ = c2 * 2.0 / (k + 1.0) - (u2 - W) * (k - 1.0) / (k + 1.0);
            Roe = ro2 * Foam::pow(c_ / c2, 2.0 / (k - 1.0));
            Ue = W - c_;
            Pe = p2 * Foam::pow(c_ / c2, 2.0 * k / (k - 1.0));
        }


    }

    if((Pe < 0.0) || (Roe < 0))
    {
        Info << "P Ro Error" << endl;
        Info << "P = " << Pe << endl;
        Info << "Ro = " << Roe << endl;
    }

    Utotal = Ue * normalVector + Vr;

    scalar convectionSpeed = Ue;
    scalar rhoState = Roe;
    vector rhoUState = Roe * Utotal;
    scalar rhoEState = Pe / (k - 1.0) + Roe * 0.5 * magSqr(Utotal);

    rhoFlux = (convectionSpeed * rhoState) * magSf;
    rhoUFlux = (convectionSpeed*rhoUState+Pe*normalVector)*magSf;
    rhoEFlux = (convectionSpeed * (rhoEState + Pe)) * magSf;
}

void Foam::kgodunovFlux::f
(
    scalar& P, scalar& p, scalar& ro, scalar& k, scalar& f, scalar& f1
) const
{
            scalar pi; 
            scalar x; 
            scalar b; 
            scalar c;
            pi = P / p;
            x = (k - 1.0) / (2.0 * k);
            b = Foam::sqrt(((k + 1.0) / (2.0 * k)) * pi + x);
            c = Foam::sqrt(k * p / ro);

            if (P<=0 || p<=0)
            {
                Info<<"P or p <= 0"<<endl;
                Info<<"P = "<<P<<endl;
                Info<<"p = "<<p<<endl;
            }
            if (pi < 0.0) 
            {
                pi *= -1.0;
                Info<<"pi <= 0"<<endl;
            }

            if (P >= p)
            {
                f = (P - p) / (ro * c * b);
                f1 = ((k + 1) * pi + 3.0 * k - 1.0) / (4.0 * ro * c * k * Foam::pow(b, 3.0));
            }
            else
            {
                f = 2.0 * c * (Foam::pow( pi, x) - 1.0) / (k - 1.0);
                f1 = c * Foam::pow(pi, x) / (k * P);
            }
}
