#pragma once
#include"wellbore1.h"
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>
#include<stdio.h>
#include<math.h>

using namespace std;

int rowA = 0;
int colA = 0;

int main(){

	////////////////////////////////////// read the data file of the well trajectoy/////////////////////////////////

	//string lineA;
	//float x;
	//string filename;
	//float arrayA[10000][3] = { { 0 } };   // define the array of the data.
	//ifstream fileIN;

	// Intro
	//cout << "input the file of well log data: measure depth, inclination angle and azimuth" << endl;
	//fileIN.open("1.txt");

	// Error Check
	//if (fileIN.fail())
	//{
	//	cout << "this file cannot be open or access" << endl;
	//}

	// reading the well log data
	//cout << "\n" << endl;
	//while (fileIN.good()){
	//	while (getline(fileIN, lineA)){
	//		istringstream streamA(lineA);
	//		colA = 0;
	//		while (streamA >> x){
	//			arrayA[rowA][colA] = x;
	//			colA++;
	//		}
	//		rowA++;
	//		}
	//	}

	// display data
	//cout << "# of row ------>" << rowA << endl;
	//cout << "# of colums----->" << colA << endl;
	//cout << " " << endl;

	//for (int i = 0; i < rowA; i++){
	//	for (int j = 0; j < colA; j++)
	//	{
	//		cout << left << setw(6) << arrayA[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	//fileIN.close();

	//double count;
	//count = rowA;
	//cout << "the number of the data points is " << count << endl;

	/////////////////////////////// input the parameter of the wellbore and formation data/////////////////////////////

	// caution!! this is field unit.

	double D_po = 5.875;           // pipe outer diameter
	//double D_pi = 1;           // pipe inner diameter                                 (unknown)
	double length_h_s1 = 6200;   // first section of the hole
	double length_h_s2 = 8365;   // second section of the hole
	double lenght_h_s3 = 1287 - 695.8;   // third section of the hole

	double TVD = 12978.54;

	double D_h_s1 = 9.56;    // first section inner diameter of the hole
	double D_h_s2 = 9.76;   // second section inner diameter of the hole
	double D_h_s3 = 9.5;    // third section inner diameter of the hole

	double e_den = 12.6;    // ppg
	double t_y = 10.06;     // yield shear stress
	double K = 0.3285;      // flow consistency
	double m = 0.7323;      // flow index
	double g = 9.8;         // gravity acceleration

	// cutting transport parameter
	double ROP = 20;                     // Check the unit ft/hour                       (back calculate)
	double ds = 2/1000000;               // Particle diameter of the      cm^2
	double h = 1/100;                    // The thickness of the cutting  1cm
	double e_cutting = 2.8;             // The density of the cutting

	// the case, you don't need to calculate the cutting!!!

	double cutting_concentration = 0;   // cuting concentration calculation model .   Feifei zhang or API 13 RD
	// Change the unit to SI unit

	// tool joint information
	double D_to = 7;                    // assume 6 inch
	double L_tooljoint = 36;            // assume 24 cm

	 D_to = D_to * 25.4 / 1000;
	 L_tooljoint = L_tooljoint / 100;         // change to meter

	//  flow rate need to calculate!

	double Q = 697.23;     // gpm            // Injection flow rate
    Q = Q * (3.785 * pow(10, -3) / 60);       // convert gpm to m^3 /s.

	 // transform inch to m     (wellbore)

	  D_h_s1 = D_h_s1 * 25.4 / 1000;    
	  D_h_s2 = D_h_s2 * 25.4 / 1000;
	  D_h_s3 = D_h_s3 * 25.4 / 1000;

	  // transform inch to m    (drillpipe)

	  D_po = D_po * 25.4 / 1000;           // pipe outer diameter

	  // transform ft to m

	  length_h_s1 = length_h_s1 * 0.3048;   // first section of the hole
	  length_h_s2 = length_h_s2 * 0.3048;   // second section of the hole
	  lenght_h_s3 = lenght_h_s3 * 0.3048;   // third section of the hole
	  TVD = TVD * 0.3048;

	  double s_length = length_h_s1 + length_h_s2 + lenght_h_s3;
	  std::cout << " the whole length of the sensor" << s_length << endl;

	 //transform lbm/gal to kg/cm^3

	 e_den = e_den * 119.83;        // change ppg to kg/cm^3 
	 //t_y = t_y * 4.44822/ 100 / 0.092903;
	 //K = K * 4.44822 / 100 / 0.092903;

	//double A_pip = PI*pow(D_pi, 2) / 4;     // drill string inner pipe diameter!!!

	double Tortuosity = 0;          // take the tortuosity into consideration..
	double Roughness_wellbore = 0;  // take roughness of wellbore into consideration..

	double Length = 0;

	//////////////////////////////////////// frictional pressure drop /////////////////////////////////////////////////

	// define an vector for storing the frictional pressure data
	vector<int> frictional_pressure_drop;

	// inner annulus pressure drop for the well
	vector<double> frictional_P_drop = { {0} };
	double frictional_drop = 0;

	// need define a step for the calculation.

	double step_length = s_length / 4000;

	for (float j = 0; j < 4000; j++)
	{
		// define the step length for each iteration
		// j is the step number

		double ToolJoint = 0;            // take tooljoint effect into consideration..

		double i = step_length * j;

		// add a check if there will be a tooljoint in the calculation.
		double drillpipe_number = 0;
		drillpipe_number  = fmod(i, 9);
		cout << drillpipe_number << "****************" << endl;

		if (drillpipe_number < 1)   // there has a tool joint in it
		{
			ToolJoint = 0;
		}

		cout << "Tooljoint Effect of the annulus " << ToolJoint << endl;

		if (i < length_h_s1)
		{
			double stdoff = 0;    // vertical well, we assume the drill pipe in concentric. 
			double D_h = D_h_s1;  // the inner diameter of the casing.
			frictional_drop = wellboreAnnulus_pressure_drop(step_length, D_po, D_to, L_tooljoint, stdoff, m, D_h, e_den, t_y, K, Q, ToolJoint, Roughness_wellbore, Tortuosity);   // pressure drop function.
		}

		if (i> length_h_s1 && i < length_h_s2)
		{
			double stdoff = 0.7;    // input the stdoff of the drill string
			double D_h = D_h_s2;     // the inner diameter of the casing
			frictional_drop = wellboreAnnulus_pressure_drop(step_length, D_po, D_to, L_tooljoint, stdoff, m, D_h, e_den, t_y, K, Q, ToolJoint, Roughness_wellbore, Tortuosity);   // pressure drop function.
		}

		if (i > length_h_s2)
		{
			double stdoff = 0.7;    // input the stdoff of the drill string
			double D_h = D_h_s3;
			frictional_drop = wellboreAnnulus_pressure_drop(step_length, D_po, D_to, L_tooljoint, stdoff, m, D_h, e_den, t_y, K, Q, ToolJoint, Roughness_wellbore, Tortuosity);   // pressure drop function.
		}

		frictional_P_drop.push_back(frictional_drop);
		frictional_P_drop[j+1] = frictional_P_drop[j] + frictional_drop;
		std::cout << " the frictional pressure drop " << frictional_P_drop[j+1] << endl;

	}
	std::system("pause");

	vector<int> gravity_pressure_drop = { {0} };
	double step = 4000;
	double step_s = TVD / step;

	double pressure_drop_gravity = 0;

	cout << "gravity pressure  " << e_den * g * TVD << endl;
	system("pause");

	for (int j = 0; j < 4000; j++)
	{
		pressure_drop_gravity = e_den * g * step_s;    // MD is the step length of the calculation
		gravity_pressure_drop.push_back(pressure_drop_gravity);
		gravity_pressure_drop[j + 1] = gravity_pressure_drop[j] + pressure_drop_gravity;
		std::cout << "  gravity pressure drop of the well , " << gravity_pressure_drop[j + 1] << endl;
	}

	std::system("pause");

	//////////////////////////////////////// gravity pressure drop //////////////////////////////////////////////////

	// caution calculate the ROP and wellbore cutting concentration!! /////  

	//////////////////////////////////////// cutting transport ratio ////////////////////////////////////////

	// calculate the cutting concentration in the wellbore
	double C_cutting = 0;
	double e_mix = 0;
	if (cutting_concentration)
	{
		if (cutting_concentration == 1)
		{
			// calculate the cutting concentration according to (Feifei Zhang. 2014)
			double A_w = PI * pow(D_h_s3, 2) / 4;
			double A_annulus = PI* (pow(D_h_s3, 2) - pow(D_po, 2)) / 4;
			// 1. superficial cutting concentration in annuli
			double v_sc = ROP * A_w / A_annulus;
			// 2. superficial liquid velocity in annuli
			double v_sl = Q / A_annulus;
			// 3. cutting feed concentration
			double CF = v_sc / (v_sc + v_sl);
			// 4. calculate the Reynold's number of the annuli (according to Ozbayoglu, 2002)
			double d_hyd = D_h_s3 - D_po;
			double d_eq = sqrt(pow(D_h_s3, 2) - pow(D_po, 2));
			double d_dim = d_hyd / d_eq;
			double a = 1.012 - 0.25* pow(d_dim, 1.563);
			double b = 0.51 - 0.242* pow(d_dim, 1.366);
			double K_prim = K * (a + b / m);
			double Re_general = pow(8, (1 - m))*e_den*pow(v_sl, (2 - m))*pow(d_hyd, m) / K_prim;
			// 5. calculate the particle slip velocity in a vertical well
			double v_slip = 0;
			if (Re_general > 1 && Re_general < 800)
			{
				double w = PI*(D_h_s3 + D_po) / 2;
				double h = (D_h_s3 - D_po) / 2;
				double A_anu = w * h;
				double v_av = Q / A_anu;
				double r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_h_s3 - D_po);     // shear rate      there is a problem!!!!  should be the generalized flow index in the formula.
				double t_w1 = t_y + K * pow(r_s, m);                                  // initial t_w1
				double e = 1000;
				do
				{
					double x = t_y / t_w1;
					double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
					double t_w2 = t_y + K * pow((r_s / Ca), m);
					e = t_w2 - t_w1;
					t_w1 = t_w2;
				} while (e > 0.00000001);

				double x = t_y / t_w1;          // the eventually data from calculation.
				double Ca = (1 - x / (1 + m) - m* pow(x, 2) / (1 + m)) * 2 / w;
				double N = Ca*m / (1 + 2 * m*(1 - Ca));
				r_s = ((1 + 2 * N) / 3 / N)*(12 * v_av) / (D_h_s3 - D_po);
				double v_slip = 0.2 * pow(9.8 *(e_cutting - e_den) / e_den, 0.72)* pow(ds, 1.18) / pow(t_w1 / r_s / e_den, 0.45);
			}
			else
			{
				double v_slip = 1.74 * pow(9.8 *(e_cutting - e_den) / e_den, 0.5)* pow(ds, 0.5);
			}
			// 6. calculate the cutting concentration of the wellbore
			double v_mix = v_sl + v_sc;
			double v_mean = (v_mix - v_slip) / 2;
			double C_cuting = v_mean + pow((pow(v_mean, 2) + v_mix*CF / v_slip), 2);
		}

		// calculate the mix flow density (According to API RD 13)
		//????????????????? transfer unit ???????????????????????????
		if (cutting_concentration == 2)
		{
			// 1. according to the correlation of Sooh
			double e_eff = h *(e_cutting - e_den) / e_den;
			double v_slip = 2.19 * pow(e_eff, 0.5);

			// 2. particle Reynolds number
			double t_s = 7.9 * pow(h*(8.345*e_cutting - e_cutting), 0.5);

			// 3. the apparent viscosity is calculated by this shear rate
			double w = PI*(D_h_s3 + D_po) / 2;
			double h = (D_h_s3 - D_po) / 2;
			double A_anu = w * h;
			double v_av = Q / A_anu;
			double r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_h_s3 - D_po);     // shear rate      there is a problem!!!!
			double t_w1 = t_y + K * pow(r_s, m);                                  // initial t_w1
			double e = 1000;
			do
			{
				double x = t_y / t_w1;
				double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
				double t_w2 = t_y + K * pow((r_s / Ca), m);
				e = t_w2 - t_w1;
				t_w1 = t_w2;
			} while (e > 0.00000001);
			double u_app = 479 * t_w1 / r_s;

			// 4. calculate the Reynolds number 
			double N_rep = 928 * e_den * v_slip;
			if (N_rep < 100)
			{
				v_slip = 0.0203*t_s * pow(ds* r_s / pow(e_den, 1 / t_w1), 0.5);
			}
			double A_annulus = PI* (pow(D_h_s3, 2) - pow(D_po, 2)) / 4;
			double v_ann = Q / A_annulus;
			double v_up = v_ann - v_slip;
			double Rt = v_up / v_ann;        // cutting transport ratio
			C_cutting = pow(D_h_s3, 2)* ROP / 448 / Q / Rt;
		}
	}

	e_mix = e_den * (1 - C_cutting) + e_cutting * C_cutting;	
	   
	// define an array for storing the gravitional pressure data.


	 //  vector<double> gravity_pressure_drop;
	 // // consider the well trajectory!     ( if we have the TVD data, we can calculate it directly)
	 // // cause we only need the final point!

	 //  for (float i = 0; i < 4619.61; i++)
	 //  {
	 //  double MD = i;
	 //  double pressure_drop_gravity = e_mix * g * MD;    // MD is the step length of the calculation
	 //  std::cout << i <<"  gravity pressure drop of the well ," << pressure_drop_gravity << endl;
	 //  }
	 //std::system("pause");
	   
	///////////////////////////////// pressure profile along the wellbore ////////////////////////////////////

	///////////////////////////////////// output the data file/////////////////////////////////////////////////////

	//	ofstream outfile;
	//	outfile.open("result the calculation");
	//	outfile << "Annulus Pressure drop" << Annulus_Pressure_Profile[i] << endl;
	//	outfile.close();
	//	return 0;
	//}
}