#ifndef ODE_DYN_H
#define ODE_DYN_H

// Include Files
#include "config_outer.h"

#define G_QUAD 9.8
#define M_QUAD 0.028

#define Kp_Z   2.0
#define Kp_VZ  25.0
#define Ki_VZ  15.0

#define Kp_rr  250.0
#define Kp_pr  250.0
#define Kp_yr  120.0

#define Ki_rr  500.0
#define Ki_pr  500.0
#define Ki_yr  16.7

#define Ip_qr	-0.760696995059653
#define Iq_pr	0.761902331719982
#define Ir_pq	-0.002866960484664
#define Im_xx	6.034380278196999e4
#define Im_yy	6.003985926176670e4
#define Im_zz	3.417442050093412e4

// define here  your ODE system yp = \dot y = f(y)
class OdeFunc {
public:
	real params[N_PARAMS];

public:
	template <class C>
	void operator()(C yp[N_STATE], C y[N_STATE]) const {
		/* yp[0] = 1 - (params[1]+1)*y[0] + params[0]*y[0]*y[0]*y[1];
		yp[1] = params[1]*y[0] - params[0]*y[0]*y[0]*y[1];*/
		C cosRoll = cos(y[0]*M_PI/180.0);
		C sinRoll = sin(y[0]*M_PI/180.0);

		C cosPitch = cos(y[1]*M_PI/180.0);
		C sinPitch = sin(y[1]*M_PI/180.0);

		C velZ_sp = Kp_Z * (params[3] - y[12]);

		C thrust = 1000.0*(Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[13]) + 36000;

		// integrale of error in p , q and r
		yp[6] = (params[0] - y[3]);//err_p;
		yp[7] = (params[1] - y[4]);//err_q;
		yp[8] = (params[2] - y[5]);//err_r;

		C cmd_r = y[6]*Ki_rr + yp[6]*Kp_rr;
		C cmd_p = y[7]*Ki_pr + yp[7]*Kp_pr;
		C cmd_y = y[8]*Ki_yr + yp[8]*Kp_yr;
		//std:cout << getAAF(cmd_p).convert_int() << std::endl;

		C temp1 = 2.7783e-12*thrust + 2.5956e-08;

		C Mx = temp1*cmd_r - 2.7783e-12*cmd_p*cmd_y;
		C My = -2.7783e-12 * cmd_r*cmd_y + temp1*cmd_p;
		C Mz = -2.5409e-13*cmd_r*cmd_p + 61.4875*temp1*cmd_y;
		C F  = 2.1354e-11 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 1.5960e-06*thrust + 0.0075;

		// Roll , pitch , yaw derivatives
		yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tan(y[1]*M_PI/180.0);
		yp[1] = y[4]*cosRoll - y[5]*sinRoll;
		yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;

		// p , q and r derivatives
		yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
		yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
		yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;

		// derivatives of body speed u , v and w 
		yp[9] = y[5]*y[10] - y[4]*y[11] + G_QUAD*sinPitch;
		yp[10] = y[3]*y[11] - y[5]*y[9] - G_QUAD*cosPitch*sinRoll;
		yp[11] = y[4]*y[9] - y[3]*y[10] + F/M_QUAD - G_QUAD*cosPitch*cosRoll;

		//Z derivative coordinate
		yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
		// Z integrale for thrust setpoint calculation
		yp[13] = velZ_sp - y[11];
	}
};

#endif // ODE_DYN_H
