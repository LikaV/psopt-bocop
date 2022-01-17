//////////////////////////////////////////////////////////////////////////
////////////////             reorientation.cxx        ////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Reorientation of an asymmetric rigid body  ///////////////
//////// Last modified: 03 July 2010                      ////////////////
//////// Reference:     Fleming et al                     ////////////////
////////                                                  ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2010        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"



//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////


adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   return 0.0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    adouble u_1 		   = controls[ 0 ];
    adouble u_2          = controls[ 1 ];
    adouble u_3          = controls[ 2 ];
    adouble q4          = controls[ 3 ];

    adouble q1 		   = states[ 0 ];
    adouble q2 		   = states[ 1 ];
    adouble q3 		   = states[ 2 ];
    adouble omega1 	   = states[ 3 ];
    adouble omega2      = states[ 4 ];
    adouble omega3      = states[ 5 ];
	adouble omega_r1	= states[ 6 ];
    	adouble omega_r2	= states[ 7 ];
    	adouble omega_r3      	= states[ 8 ];
    	

        double Ix = 5621.0;
        double Iy = 4547.0;
        double Iz = 2364.0;

	double J_1 = Ix;
    	double J_2 = Iy;
   	double J_3 = Iz;
	double J_r1 = 115.0;
    	double J_r2 = 94;
    	double J_r3 = 53;

        adouble dq1     = 0.5*( omega1*q4 - omega2*q3 + omega3*q2 );
        adouble dq2     = 0.5*( omega1*q3 + omega2*q4 - omega3*q1 );
        adouble dq3     = 0.5*(-omega1*q2 + omega2*q1 + omega3*q4 );
    
	adouble domega1 = 1/J_1*((J_2*omega2+J_r2*omega_r2)*omega3-(J_3*omega3+J_r3*omega_r3)*omega2-J_r1*u_1);
    	adouble domega2 = 1/J_2*((J_3*omega3+J_r3*omega_r3)*omega1-(J_1*omega1+J_r1*omega_r1)*omega3-J_r2*u_2);
    	adouble domega3 = 1/J_3*((J_1*omega1+J_r1*omega_r1)*omega2-(J_2*omega2+J_r2*omega_r2)*omega1-J_r3*u_3);
	adouble domega_r1 = u_1;
    	adouble domega_r2 = u_2;
    	adouble domega_r3 = u_3;




    //adouble domega1 = u1/Ix - ((Iz-Iy)/Ix)*omega2*omega3;
    //adouble domega2 = u2/Iy - ((Ix-Iz)/Iy)*omega1*omega3;
    //adouble domega3 = u3/Iz - ((Iy-Ix)/Iz)*omega1*omega2;




    derivatives[ 0 ] =   dq1;
    derivatives[ 1 ] =   dq2;
    derivatives[ 2 ] =   dq3;
    derivatives[ 3 ] =   domega1;
    derivatives[ 4 ] =   domega2;
    derivatives[ 5 ] =   domega3;
    derivatives[ 6 ] =   domega_r1;
    derivatives[ 7 ] =   domega_r2;
    derivatives[ 8 ] =   domega_r3;


    path[ 0 ] = q1*q1 + q2*q2 + q3*q3 + q4*q4 - 1.0;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

    adouble q1i 	      = initial_states[ 0 ];
    adouble q2i 	      = initial_states[ 1 ];
    adouble q3i 	      = initial_states[ 2 ];
    adouble omega1i 	   = initial_states[ 3 ];
    adouble omega2i     = initial_states[ 4 ];
    adouble omega3i     = initial_states[ 5 ];
    adouble omegar1i 	   = initial_states[ 6 ];
    adouble omegar2i     = initial_states[ 7 ];
    adouble omegar3i     = initial_states[ 8 ];

    adouble initial_controls[4], final_controls[4], q4i, q4f;

    get_initial_controls( initial_controls, xad, iphase, workspace );

    get_final_controls(   final_controls  , xad, iphase, workspace );

    q4i = initial_controls[ 3 ];

    q4f = final_controls[   3 ];

    adouble q1f 	      = final_states[ 0 ];
    adouble q2f 	      = final_states[ 1 ];
    adouble q3f 	      = final_states[ 2 ];
    adouble omega1f 	   = final_states[ 3 ];
    adouble omega2f     = final_states[ 4 ];
    adouble omega3f     = final_states[ 5 ];
    adouble omegar1f 	   = initial_states[ 6 ];
    adouble omegar2f     = initial_states[ 7 ];
    adouble omegar3f     = initial_states[ 8 ];


    e[ 0  ] 	= q1i;
    e[ 1  ]    = q2i;
    e[ 2  ]    = q3i;
    e[ 3  ]   	= q4i;
    e[ 4  ]   	= omega1i;
    e[ 5  ]    = omega2i;
    e[ 6  ]  	= omega3i;
    e[ 7  ]   	= omegar1i;
    e[8  ]    = omegar2i;
    e[ 9  ]  	= omegar3i;

    e[ 10 ] 	= q1f;
    e[ 11  ]    = q2f;
    e[ 12]     = q3f;
    e[ 13 ]    = q4f;
    e[ 14 ]   	= omega1f;
    e[ 15 ]    = omega2f;
    e[ 16 ]  	= omega3f;
    e[ 17 ]   	= omegar1f;
    e[ 18 ]    = omegar2f;
    e[ 19 ]  	= omegar3f;



}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{

   // Single phase problem

}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Spacecraft reorientation";
    problem.outfilename                 = "reorientation.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   					= 1;
    problem.nlinkages               = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 9;
    problem.phases(1).ncontrols 		= 4;
    problem.phases(1).nevents   		= 20;
    problem.phases(1).npath     		= 1;
    problem.phases(1).nodes         << 40;


    psopt_level2_setup(problem, algorithm);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states << -1.0,    -1.0,    -1.0,   -0.5,    -0.5,  -0.5, -4.5,    -4.5,  -4.5;
    problem.phases(1).bounds.upper.states <<  1.0,     1.0,     1.0,    0.5,     0.5,   0.5,  4.5,     4.5,   4.5;


    problem.phases(1).bounds.lower.controls << -3.0,  -3.0, -3.0,  -1.0;

    problem.phases(1).bounds.upper.controls <<  3.0,   3.0,  3.0,   1.0;


    problem.phases(1).bounds.lower.path(0) = 0.000;

    problem.phases(1).bounds.upper.path(0) = 0.000;

    MatrixXd q0(4,1), qf(4,1), omega0(3,1), omegaf(3,1);//, omegar0(3,1), omegarf(3,1);

    double phi = 150.0*(pi/180.0);

    q0(0) = 0.0;          q0(1) = 0.0; q0(2) = 0.0; q0(3) = 1.0;
    qf(0) = sin(phi/2.0); qf(1) = 0.0; qf(2) = 0.0; qf(3) = cos(phi/2.0);


    omega0 = zeros(1,3); omegaf = zeros(1,3);
     

    problem.phases(1).bounds.lower.events(0) = q0(0); 
    problem.phases(1).bounds.lower.events(1) = q0(1);
    problem.phases(1).bounds.lower.events(2) = q0(2);
    problem.phases(1).bounds.lower.events(3) = q0(3);     
    problem.phases(1).bounds.lower.events(4) = omega0(0);
    problem.phases(1).bounds.lower.events(5) = omega0(1);
    problem.phases(1).bounds.lower.events(6) = omega0(2); 

    problem.phases(1).bounds.lower.events(7) = omega0(0);
    problem.phases(1).bounds.lower.events(8) = omega0(1);
    problem.phases(1).bounds.lower.events(9) = omega0(2);  
    
    
    problem.phases(1).bounds.lower.events(10)  = qf(0); 
    problem.phases(1).bounds.lower.events(11)  = qf(1);
    problem.phases(1).bounds.lower.events(12)  = qf(2);  
    problem.phases(1).bounds.lower.events(13) = qf(3);  
    problem.phases(1).bounds.lower.events(14) = omegaf(0);
    problem.phases(1).bounds.lower.events(15) = omegaf(1);
    problem.phases(1).bounds.lower.events(16) = omegaf(2);

    problem.phases(1).bounds.lower.events(17) = omegaf(0);
    problem.phases(1).bounds.lower.events(18) = omegaf(1);
    problem.phases(1).bounds.lower.events(19) = omegaf(2);            
    
   
    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;
    

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.01;
    problem.phases(1).bounds.upper.EndTime      = 50.0;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 		= &integrand_cost;
    problem.endpoint_cost 			= &endpoint_cost;
    problem.dae             		= &dae;
    problem.events 					= &events;
    problem.linkages					= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			= problem.phases(1).nodes(0);
    int ncontrols          = problem.phases(1).ncontrols;
    int nstates            = problem.phases(1).nstates;

    MatrixXd state_guess    =  zeros(nstates,nnodes);
    MatrixXd control_guess  =  2.0*ones(ncontrols,nnodes);
    MatrixXd time_guess     =  linspace(0.0,45,nnodes);


    state_guess              << linspace( q0(0), qf(0), nnodes ),
                                linspace( q0(1), qf(1), nnodes ),
                                linspace( q0(2), qf(2), nnodes ),
                                linspace( omega0(0), omegaf(0), nnodes ),
                                linspace( omega0(1), omegaf(1), nnodes ),
                                linspace( omega0(2), omegaf(2), nnodes ),
				linspace( omega0(0), omegaf(0), nnodes ),
                                linspace( omega0(1), omegaf(1), nnodes ),
                                linspace( omega0(2), omegaf(2), nnodes );
                                
    control_guess.row(3) = linspace( q0(3), qf(3), nnodes );

    problem.phases(1).guess.states   = state_guess;
    problem.phases(1).guess.controls = control_guess;
    problem.phases(1).guess.time     = time_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    algorithm.ode_tolerance               = 1.e-5;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd states, controls, t, q1, q2, q3, q4, omega, u, q,omegar;

    states      = solution.get_states_in_phase(1);
    controls    = solution.get_controls_in_phase(1);
    t           = solution.get_time_in_phase(1);

    q1 = states.row(0);
    q2 = states.row(1);
    q3 = states.row(2);
    q4 = controls.row(3);
    
    q.resize(4,length(t));

    q << q1 , 
         q2 , 
         q3 , 
         q4;

    omega.resize(3, length(t)); 
    
    omega << states.row(3),
             states.row(4),
             states.row(5);

    omegar.resize(3, length(t)); 

    omegar << states.row(6),
             states.row(7),
             states.row(8);
    
    u.resize(3,length(t));
     
    u << controls.row(0),
         controls.row(1),
         controls.row(2);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(states, "states.dat");
    Save(controls, "controls.dat");
    Save(t,"t.dat");



////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////



    multiplot(t,q,problem.name+": quarternion","time (s)", "q1 q2 q3 q4", "q1 q2 q3 q4", 2, 2);

    multiplot(t,u, problem.name+": controls", "time (s)", "u1 u2 u3", "u1 u2 u3");

    multiplot(t,omega, problem.name+": omega",  "time(s)","omega1 omega2 omega3","omega1 omega2 omega3",3,1);
    multiplot(t,omegar, problem.name+": omegar",  "time(s)","omegar1 omegar2 omegar3","omegar1 omegar2 omegar3",3,1);

    multiplot(t,q,problem.name+": quarternion","time (s)", "q1 q2 q3 q4", "q1 q2 q3 q4", 2, 2,
                "pdf", "reorientation_q.pdf");

    multiplot(t,u, problem.name+": controls", "time (s)", "u1 u2 u3", "u1 u2 u3", 3, 1,
               "pdf", "reorientation_u.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
