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

#include "agodunovFlux.H"

//#include <stdio.h>
//#include <time.h>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::agodunovFlux::evaluateFlux
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
//вектор нормали
const vector normalVector = Sf / magSf;
//давление
scalar p_l = pLeft;
scalar p_p = pRight;
//плотность
scalar ro_l = pLeft / (RLeft * TLeft);
scalar ro_p = pRight / (RRight * TRight);
//скорость звука
scalar c_l, c_p;
//показатель адиабаты
scalar kappa_p = (CvRight + RRight) / CvRight;
scalar kappa_l = (CvLeft + RLeft) / CvLeft;

/*
 * Разбиваем вектор скорости на компоненты - нормальную и тангенциальную
 * В случае 2-мера Utau - перпендикулярная компонента скорости к вектору нормали грани ячейки 
 * В случае 3-мера Utau - сумма 2 компонент вектора скорости, ортогональных вектору нормали поверхности ячейки
*/
vector Unorm1, Unorm2;
vector Utau1, Utau2;
//Скорости УВ, ВР, Вакуума
scalar U_ud, U_raz, U_vak;
//разделение скорости
Unorm1 = (ULeft & normalVector) * normalVector;
Unorm2 = (URight & normalVector) * normalVector;
Utau1 = ULeft - Unorm1;
Utau2 = URight - Unorm2;

vector U_l = Unorm1;
vector U_p = Unorm2;

vector V_l = Utau1;
vector V_p = Utau2;

//массовые скорости
scalar a_l, a_p;
vector D_l, D_p, D_le, D_pe, c_le, c_pe;
scalar R_l, R_p;
/*
 * P - БОЛЬШОЕ ДАВЛЕНИЕ из Итерационного процесса (метод Ньютона)
 * p0 - минимальное давление ВАКУУМА
 * K1 - переменная скорости для конфигурации ВР слева
 * K2 - переменная скорости для конфигурации ВР справа
 * K4 - i-ое приращение БОЛЬШОГО ДАВЛЕНИЯ в Итерационном процессе
 */
scalar P, p0, K1, K2, K4;
//в методе Ньютона К4 - изменение давления на i-ом шаге
K4 = 1000000.;
scalar f_p, f_l, f1_p, f1_l;

//Большие величины
scalar Pr, Rr;
vector U, Ur, Vr, Utotal;
/*
 * Показатель адибаты, используемый для вычисления консервативных переменных, 
 * определяется из сравнения тангенциальной компоненты скорости с нулём
*/	 
scalar k;

//поправки на вакуум и скорости движения грани ячейки
scalar p0_l = 0.;
scalar p0_p = 0.;
scalar W = 0.;
/*
 * Если давление слева меньше давления справа, то не меняем
*/
if( p_l <= p_p) {

	c_l = sqrt(kappa_l * (p_l + p0_l) / ro_l);
	c_p = sqrt(kappa_p * (p_p + p0_p) / ro_p);
	U_ud = (p_p - p_l) / sqrt(ro_l*( (kappa_l + 1.0) * (p_p + p0_l) / 2.0 + (kappa_l - 1.0) * (p_l+p0_l)/  2.0) );
	U_raz = -2 * c_p * (1.0 - pow( ( (p_l + p0_p) / (p_p + p0_p) ),( ( kappa_p - 1.0) / (2.0 * kappa_p) ) ) ) / (kappa_p - 1.0);    

	//Выбираем минимальное нулевое давление 
	if (p0_p <= p0_l) {
		U_vak = -(2.0 * c_l / (kappa_l - 1.0)) * (1.0 - pow( ( (p0_l - p0_p) / (p_l + p0_l) ),( (kappa_l - 1.0) / (2.0 * kappa_l) ) )) - 2.0 * c_p / (kappa_p - 1.0);
	} else {
		U_vak = -(2.0 * c_p / (kappa_p - 1.0)) * (1.0 - pow( ( (p0_p - p0_l) / (p_p + p0_p) ),( (kappa_p - 1.0) / (2.0 * kappa_p) ) )) - 2.0 * c_l / (kappa_l - 1.0); 
	}

	if ( ((U_l - U_p) & normalVector) > U_raz) {//Конфигурация  > ВР

	if (fabs( ((U_l - U_p) & normalVector) - U_ud) < 1.e-3) {
		P = p_p;// Начальное приближение
	} else {
		P = (p_l * ro_p * c_p + p_p * ro_l * c_l + ((U_l - U_p) & normalVector) * ro_l * c_l * ro_p * c_p) / (ro_p * c_p + ro_l * c_l);// Начальное приближение УВ
		if ((P - 0.) < 1.e-3) P = p_p;// Начальное приближение если 
		int count = 0;
		
		while (fabs( K4 / P) > 1.e-6) {
			count++;
			if (((U_l-U_p) & normalVector) > U_ud) {//две ударные волны
				f_l = (P - p_l) / (ro_l * c_l * sqrt( (kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l) )); //ударная
				f1_l = ((kappa_l + 1.0) * (P + p0_l) / (p_l + p0_l) + (3.0 * kappa_l - 1.0)) / (4.0 * kappa_l * ro_l * c_l * pow( ((kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l)), (3/2)) ); //ударная

				f_p = (P - p_p) / (ro_p * c_p * sqrt( (kappa_p + 1.0) * (P + p0_p) / (2.0 * kappa_p * (p_p + p0_p)) + (kappa_p - 1.0) / (2.0 * kappa_p) )); //ударная
				f1_p = ((kappa_p + 1.0) * (P + p0_p) / (p_p + p0_p) + (3.0 * kappa_p - 1.0)) / (4.0 * kappa_p * ro_p * c_p * pow( ((kappa_p + 1.0) * (P + p0_p) / (2.0 * kappa_p * (p_p + p0_p)) + (kappa_p - 1.0) / (2.0 * kappa_p)), (3/2)) ); //ударная

			} else {//левая ударная волна и правая волна разряжения

				f_l = (P - p_l) / (ro_l * c_l * sqrt( (kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l) )); //ударная
				f1_l = ((kappa_l + 1.0) * (P + p0_l) / (p_l + p0_l) + (3.0 * kappa_l - 1.0)) / (4.0 * kappa_l * ro_l * c_l * pow( ((kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l)), (3/2)) ); //ударная

				f_p = 2.0 * c_p * ( pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) - 1.0) / (kappa_p - 1.0); //правая волна разрежения
				f1_p = c_p * pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) / (kappa_p * (P + p0_p)); //правая волна разрежения
			}
			
			if (count > 1000) {
				Info << "Превышено количество итераций 1. Р = " << P <<  " K4 = " << K4 << endl;
			}
			
			K4 = f_l + f_p - ((U_l - U_p) & normalVector);
			P = P - K4 / (f1_l + f1_p); //итеррационный процесс
			if ((P - 0.) < 1.e-7) break;
		}
	}
	} else if (((U_l - U_p) & normalVector) > U_vak) {//Конфигурация > Вакуума

		if(fabs(((U_l - U_p) & normalVector) - U_raz) < 1.e-3) {
			P = p_l;// Начальное приближение
		} else {
			P = (p_l * ro_p * c_p + p_p * ro_l * c_l + ((U_l - U_p) & normalVector) * ro_l * c_l * ro_p * c_p) / (ro_p * c_p + ro_l * c_l);// Начальное приближение ВР      
		if ((P - 0.) < 1.e-3) P = p_l;
		int count = 0;
		
		while (fabs(K4 / P) > 1.e-6) {
			count++;
			
			f_l = 2.0 *c_l * (pow( ((P + p0_l) / (p_p + p0_p)),((kappa_l - 1.0) / (2.0 * kappa_l)) ) - 1.0) / (kappa_l - 1.0); //левая волна разрежения
			f1_l = c_l * pow( ((P + p0_l) / (p_p + p0_l)),((kappa_l - 1.0) / (2.0 * kappa_l)) ) / (kappa_l * (P + p0_l)); //левая волна разрежения
			
			f_p = 2 * c_p * (pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) - 1.0) / (kappa_p - 1.0); //правая волна разрежения
			f1_p = c_p * pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) / (kappa_p * (P + p0_p)); //правая волна разрежения

			if (count > 1000) {
			Info << "Превышено количество итераций 2. Р = " << P <<  " K4 = " << K4 << endl;
			}

			K4 = f_l + f_p - ((U_l - U_p) & normalVector);
			P = P - K4 / (f1_l + f1_p); //итеррационный процесс

			if ((P - 0.) < 1.e-7) break;
		}

		}
	} else {//Вакуум
		if(p0_l >= p0_p) {
			p0 = p0_l;
		} else {
			p0 = p0_p;
		}

		if(fabs(((U_l - U_p) & normalVector) - U_vak) < 0.001) {
			P = 100 - p0;
		} else {
			P = 0 - p0;
		}
	}

	a_l = sqrt(ro_l * ((kappa_l + 1.0)*(P + p0_l) / 2.0 + (kappa_l - 1.0) * (p_l + p0_l) / 2.0));
	a_p = sqrt(ro_p * ((kappa_p + 1.0)*(P + p0_p) / 2.0 + (kappa_p - 1.0) * (p_p + p0_p) / 2.0));

} else {
	//Меняем ячейки местами
	vector u1_ = U_l;
	vector u2_ = U_p;
	scalar p1_ = p_l;
	scalar p2_ = p_p;
	scalar ro1_ = ro_l;
	scalar ro2_ = ro_p;
	scalar k1_ = kappa_l;
	scalar k2_ = kappa_p;
	scalar p01_ = p0_l;
	scalar p02_ = p0_p;

	U_l = - u2_;
	U_p = - u1_;
	p_l = p2_;
	p_p = p1_;
	ro_l = ro2_;
	ro_p = ro1_;
	kappa_l = k2_;
	kappa_p = k1_;
	p0_l = p02_;
	p0_p = p01_;


	c_l = sqrt(kappa_l * (p_l + p0_l) / ro_l);

	c_p = sqrt( kappa_p * (p_p + p0_p) / ro_p);

	U_ud = (p_p - p_l) / sqrt(ro_l * ((kappa_l + 1.0) * (p_p + p0_l) / 2.0 + (kappa_l - 1.0) * (p_l + p0_l) / 2.0));

	U_raz = -2.0 * c_p * (1.0 - pow( ((p_l + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) )) / (kappa_p - 1.0);

	if (p0_p <= p0_l) {
		U_vak = - (2.0 * c_l / (kappa_l - 1.0)) * (1.0 - pow( ((p0_l - p0_p) / (p_l + p0_l)),((kappa_l - 1.0) / (2.0 * kappa_l)) )) - 2.0 * c_p / (kappa_p - 1.0);
	} else {
		U_vak = - (2.0 * c_p / (kappa_p - 1.0)) * (1.0 - pow( ((p0_p - p0_l) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) )) - 2.0 * c_l / (kappa_l - 1.0); 
	}

	if (((U_l - U_p) & normalVector) > U_raz) {//Конфигурация > ВР

		if (fabs(((U_l - U_p) & normalVector) - U_ud) < 1.e-3) {
			P = p_p;
		} else {
			P = (p_l * ro_p * c_p + p_p * ro_l * c_l + ((U_l - U_p) & normalVector) * ro_l * c_l * ro_p * c_p) / (ro_p * c_p + ro_l * c_l);
			if ((P - 0.) < 1.e-3) P = p_p;
			int count = 0;
			
			while (fabs(K4 / P) > 1.e-6) {
				count++;
				if (((U_l - U_p) & normalVector) > U_ud) {//две ударные волны
					f_l = (P - p_l) / (ro_l * c_l * sqrt( (kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l) )); //ударная
					f1_l = ((kappa_l + 1.0) * (P + p0_l) / (p_l + p0_l) + (3.0 * kappa_l - 1.0)) / (4.0 * kappa_l * ro_l * c_l * pow( ((kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l)), (3/2)) ); //ударная

					f_p = (P - p_p) / (ro_p * c_p * sqrt( (kappa_p + 1.0) * (P + p0_p) / (2.0 * kappa_p * (p_p + p0_p)) + (kappa_p - 1.0) / (2.0 * kappa_p) )); //ударная
					f1_p = ((kappa_p + 1.0) * (P + p0_p) / (p_p + p0_p) + (3.0 * kappa_p - 1.0)) / (4.0 * kappa_p * ro_p * c_p * pow( ((kappa_p + 1.0) * (P + p0_p) / (2.0 * kappa_p * (p_p + p0_p)) + (kappa_p - 1.0) / (2.0 * kappa_p)), (3/2)) ); //ударная
				} else {   //левая ударная волна и правая волна разряжения
					f_l = (P - p_l) / (ro_l * c_l * sqrt( (kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l) )); //ударная
					f1_l = ((kappa_l + 1.0) * (P + p0_l) / (p_l + p0_l) + (3.0 * kappa_l - 1.0)) / (4.0 * kappa_l * ro_l * c_l * pow( ((kappa_l + 1.0) * (P + p0_l) / (2.0 * kappa_l * (p_l + p0_l)) + (kappa_l - 1.0) / (2.0 * kappa_l)), (3/2)) ); //ударная

					f_p = 2.0 * c_p * ( pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) - 1.0) / (kappa_p - 1.0); //правая волна разрежения
					f1_p = c_p * pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) / (kappa_p * (P + p0_p)); //правая волна разрежения
				}
				
				if (count > 1000) {
					Info << "Превышено количество итераций 3. Р = " << P <<  " K4 = " << K4 << endl;
				}

				K4 = f_l + f_p - ((U_l-U_p) & normalVector);
				P = P - K4 / (f1_l + f1_p); //итеррационный процесс

				if ((P - 0.) < 1.e-7) break;
			}
		}
	} else {

		if (((U_l - U_p) & normalVector) > U_vak) {//Конфигурация > Вакуума

			if(fabs(((U_l - U_p) & normalVector) - U_raz) < 1.e-3) {
				P = p_l;
			} else {
				P = (p_l * ro_p * c_p + p_p * ro_l * c_l + ((U_l - U_p) & normalVector) * ro_l * c_l * ro_p * c_p) / (ro_p * c_p + ro_l * c_l);
				if ((P - 0.) < 1.e-3) P = p_l;
				int count = 0;
				
				while (fabs(K4 / P) > 1.e-6) {
				count++;
				
				f_l = 2.0 *c_l * (pow( ((P + p0_l) / (p_p + p0_p)),((kappa_l - 1.0) / (2.0 * kappa_l)) ) - 1.0) / (kappa_l - 1.0); //левая волна разрежения
				f1_l = c_l * pow( ((P + p0_l) / (p_p + p0_l)),((kappa_l - 1.0) / (2.0 * kappa_l)) ) / (kappa_l * (P + p0_l)); //левая волна разрежения
				
				f_p = 2 * c_p * (pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) - 1.0) / (kappa_p - 1.0); //правая волна разрежения
				f1_p = c_p * pow( ((P + p0_p) / (p_p + p0_p)),((kappa_p - 1.0) / (2.0 * kappa_p)) ) / (kappa_p * (P + p0_p)); //правая волна разрежения

				if (count > 1000) {
					Info << "Превышено количество итераций 4. Р = " << P <<  " K4 = " << K4 << endl;
				}

				K4 = f_l + f_p - ((U_l - U_p) & normalVector);
				P = P - K4 / (f1_l + f1_p); //итеррационный процесс
				if ((P - 0.) < 1.e-7) break;
				}
			}
		} else {//Вакуум

			if(p0_l >= p0_p) {
				p0 = p0_l;
			} else {
				p0 = p0_p;
			}

			if(fabs(((U_l-U_p) & normalVector) - U_vak) < 0.001) {
				P = 100 - p0;
			} else {
				P = 0 - p0;
			}
		}
	}

	//После определения "БОЛЬШОГО" давления, возвращаем ячейки на место
	U_l = u1_;
	U_p = u2_;
	p_l = p1_;
	p_p = p2_;
	ro_l = ro1_;
	ro_p = ro2_;
	kappa_l = k1_;
	kappa_p = k2_;
	p0_l = p01_;
	p0_p = p02_;

	c_l = sqrt(kappa_l * (p_l + p0_l) / ro_l);
	c_p = sqrt(kappa_p * (p_p + p0_p) / ro_p);
	a_l = sqrt(ro_l * ((kappa_l + 1.0) * (P + p0_l) / 2.0 + (kappa_l - 1.0) * (p_l + p0_l) / 2.0));
	a_p = sqrt(ro_p * ((kappa_p + 1.0) * (P + p0_p) / 2.0 + (kappa_p - 1.0) * (p_p + p0_p) / 2.0));
} 
//Скорость КР
U = (a_l * U_l + a_p * U_p + (p_l - p_p) * normalVector) / (a_l + a_p);

if(P >= p_l) {
	D_l = U_l - (a_l / ro_l * normalVector);
	D_le = D_l;
	R_l = ro_l * a_l / (a_l - ro_l * ((U_l - U) & normalVector));
} else {
	c_le = c_l * normalVector + (kappa_l - 1) * (U_l - U) / 2;
	D_l = U_l - c_l * normalVector;
	D_le = U - c_le;
	R_l = kappa_l * (P + p0_l) / ((c_le & normalVector) * (c_le & normalVector));
}

if(P >= p_p) {
	D_p = U_p + (a_p / ro_p * normalVector);
	D_pe = D_p;
	R_p = ro_p * a_p / (a_p + ro_p * ((U_p - U) & normalVector));
} else {
	c_pe = c_p * normalVector - (kappa_p - 1.0) * (U_p - U) / 2;
	D_p = U_p + c_p * normalVector;
	D_pe = U + c_pe;
	R_p = kappa_p * (P + p0_p) / ((c_pe & normalVector) * (c_pe & normalVector));
}
//Область до распада слева
if( (D_l & normalVector) > W) {
	Ur = U_l;
	Pr = p_l;
	Rr = ro_l;
}
//Область до распада справа
if( (D_p & normalVector) < W) {
	Ur = U_p;
	Pr = p_p;
	Rr = ro_p;
}
//УВ_УВ
if((P >= p_l) && ( (D_l & normalVector) <= W) && ( (U & normalVector) >= W)) {
	Ur = U;
	Pr = P;
	Rr = R_l;
}
//УВ_УВ
if((P >= p_p) && ( (U & normalVector) <= W) && ( (D_p & normalVector) >= W)) {
	Ur = U;
	Pr = P;
	Rr = R_p;
}
//УВ_ВР
if((P <= p_l) && ( (D_le & normalVector) <= W) && ( (U & normalVector) >= W)) {
	Ur = U;
	Pr = P;
	Rr = R_l;
}
//ВР-ВР
if((P <= p_l) && ( (D_l & normalVector) <= W) && ( (D_le & normalVector) >= W)) {
	K1 = ((kappa_l - 1.0) * ((U_l & normalVector) - W) / (kappa_l + 1.0)) + 2.0 * c_l / (kappa_l + 1.0);
	Ur = (W + K1) * normalVector;
	Pr = (p_l + p0_l) * pow( (K1 / c_l),(2.0 * kappa_l / (kappa_l - 1.0)) ) - p0_l;
	Rr = kappa_l * (P + p0_l) / (K1 * K1);
}
//УВ_ВР
if((P <= p_p) && ( (D_pe & normalVector) >= W) && ( (U & normalVector) <= W)) {
	Ur = U;
	Pr = P;
	Rr = R_p;
}
//ВР-ВР
if((P <= p_p) && ( (D_pe & normalVector) <= W) && ( (D_p & normalVector) >= W)) {
	K2 = (2.0 * c_p / (kappa_p + 1.0)) - (kappa_p - 1.0) * ((U_p & normalVector) - W) / (kappa_p + 1.0);
	Ur = (W - K2) * normalVector;
	Pr = (p_p + p0_p) * pow( (K2 / c_p),(2.0 * kappa_p / (kappa_p - 1.0)) ) - p0_p;
	Rr = kappa_p * (P + p0_p) / (K2 * K2);
}

if(p0_l >= p0_p) {
	p0 = p0_l;
} else {
	p0 = p0_p;
}
//Область вакуума
if((P <= -p0) && ( (D_pe & normalVector) >= W) && ( (D_le & normalVector) <= W)) {
	//Info << "Вакуум" << endl;
	Pr = -p0;
	Ur = U;
	Rr = 0;
}

if ( (Ur & normalVector) > 0.) {//Невозмущенна среда слева
	Vr = V_l;
	k = kappa_l;
} else {//Невозмущенна среда справа
	Vr = V_p;
	k = kappa_p;
}

Utotal = Ur + Vr;

scalar convectionSpeed = (Utotal & normalVector);
scalar rhoState = Rr;
vector rhoUState = Rr * Utotal;
scalar rhoEState = Rr * (Pr / Rr / ( k - 1) + 0.5 * magSqr(Utotal));

rhoFlux = (convectionSpeed * rhoState) * magSf;
rhoUFlux = (convectionSpeed * rhoUState + Pr * normalVector) * magSf;
rhoEFlux = (convectionSpeed * (rhoEState + Pr)) * magSf;

}
