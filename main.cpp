#include "sparselizard.h"

#include "plasticity.h"
#include "loading.h"


using namespace sl;
// ----------------------------------------------------------------------------
// mesh making function
mesh createmesh(double w, double len,
                vector<double> thick_lays,
                int nx, int nz, vector<int> num_nodes_lays);


// field concatenating function
std::vector<field> concatenate(std::vector<std::vector<field>> input);

// cooling-down temperatures vector maker
std::vector<double> make_Tvec(double T_0, double T_f, double dT);

// force (axial disp & transverse force) vector maker
std::vector<double> make_Fvec(double d0, double dd, int steps);

// ----------------------------------------------------------------------------
// declaring pointer based vectors for storing strain,stress & ep in layers
std::vector<double>* L_strain_xx = NULL;
std::vector<double>* L_strain_yy = NULL;
std::vector<double>* L_strain_zz = NULL;
std::vector<double>* L_strain_yz = NULL;
std::vector<double>* L_strain_xz = NULL;
std::vector<double>* L_strain_xy = NULL;

std::vector<double>* L_stress_xx = NULL;
std::vector<double>* L_stress_yy = NULL;
std::vector<double>* L_stress_zz = NULL;
std::vector<double>* L_stress_yz = NULL;
std::vector<double>* L_stress_xz = NULL;
std::vector<double>* L_stress_xy = NULL;

// "main vector of pointers" for storing "pointer vectors"
std::vector<vector<double>*> strain_xx, strain_yy, strain_zz, strain_yz, strain_xz, strain_xy;
std::vector<vector<double>*> stress_xx, stress_yy, stress_zz, stress_yz, stress_xz, stress_xy;


// vectors for whole tape
std::vector<double> tape_strain_xx, tape_stress_xx;
std::vector<double> tape_strain_yy, tape_stress_yy;
std::vector<double> tape_strain_zz, tape_stress_zz;
std::vector<double> tape_strain_yz, tape_stress_yz;
std::vector<double> tape_strain_xz, tape_stress_xz;
std::vector<double> tape_strain_xy, tape_stress_xy;

// temperature vector
std::vector<double> vecT;
// disp. loading vector
std::vector<double> vecDispLoad;
// Tr. force vector
std::vector<double> vecTrForce;


// function for making "pointer vectors of layers" in "main vector of pointers"
void prepvecs(int num_lays);

// data extractor (used for strain_tot, stress_corr & ep)
void extractdata(mesh mymesh, int body, int avg_box, double vol,
            expression strain_tot, expression stress_corr,
            vector<double> thick_lays,
            double w, double len);



// save all the data in vectors into text files
void savedata(vector<string> tags);

// *****************************************************************************
// *****************************************************************************
int main(void)
{
    // timer
    wallclock myclock;

    // physical regions numbers
    int Copper = 1, Silver = 2, Hastelloy = 3, Buffer = 4, YBCO = 5; // materials
    int body = 6;           // body : all tapes
    int st = 7, en = 8;     // tape's vertical surface at either sides
    int bot = 9, top = 10;  // bottom & top surfaces of the tape
    int mid_top = 11;       // top surface region for transverse force
    int avg_box = 12;       // vol. for avg. tape
    int cent_line = 13;     // center line for avg. behavior for transverse loading
    int mid_box = 14;       // the box subjected to the trans. load


    // tape dimensions (in [m])
    double w = 5.9e-3, len = 4e-2; // width & length
    //double thCu = 2e-5, thHastelloy = 0.1e-3, thYBCO = 2e-6; // layers' thickness
    double thCu = 2e-5, thAg = 2e-6, thHastelloy = 0.1e-3, thBuffer = 1e-6, thYBCO = 2e-6; // layers' thickness

    // layers info
    int num_lays = 7; // number of layers
    //std::vector<double> thick_lays = {thCu, thHastelloy, thYBCO, thCu};
    //std::vector<string> tags = {"botcu", "sub", "YBCO", "topcu"}; // tags for layers
    std::vector<double> thick_lays = {thCu, thAg, thHastelloy, thBuffer, thYBCO, thAg, thCu};
    std::vector<string> tags = {"botcu", "botag", "sub", "buffer", "YBCO", "topag", "topcu"}; // tags for layers


    // number of nodes in mesh
    int nx       = 9;   // nodes across the width of the tape
    int nyCu     = 3;   // nodes along the height of Copper layer
    int nySub    = 4;   // nodes along the height of Hastelloy layer
    int nyYBCO   = 2;   // nodes along the height of YBCO layer
    int nyAg     = 2;   // nodes along the height of Silver layer
    int nyBuffer = 2;   // nodes along the height of Buffer layer
    int nz       = 101; // nodes along the length of the tape
    std::vector<int> num_nodes_lays = {nyCu, nyAg, nySub, nyBuffer, nyYBCO, nyAg, nyCu};


    // creating the mesh with the function
    mesh mymesh = createmesh(w, len,
                              thick_lays,
                              nx, nz, num_nodes_lays);


    // volume of the body
    double vol = expression(1).integrate(avg_box, 4);
    double vol_mid_box = expression(1).integrate(mid_box, 4);


    // making vectors for data
    prepvecs(num_lays);

    // Simulation main settings
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //
    //int study_case = 0; // elastic
    int study_case = 1; // plastic


    // (***************)
    // Thermal loading
    bool th_load = false;

    double T_room = 300.0;   // [k]
    double T_final = 77.0;   // [K]
    double dT = 1.0;

    vecT = make_Tvec(T_room, T_final, dT);   // vector of temps to cool-down
    int steps_I = vecT.size();               // steps of cooling-down

    // final delta_T for last step of cool-down
    double dT_final = vecT[steps_I-1] - vecT[steps_I-2];
    // !!! we may have one last different delta_T (e.g. 300 >> 77, 2 for all and 1 for the last)

    // delta_T as parameter (for sig_thermal)
    parameter delta_T;


    // (***************)
    // axial loading (zz)
    bool axial_disp = false;

    int steps_II = 399;     // steps of loading
    double disp0 = 1.0 * 1e-6;
    double delta_disp = 1.0 * 1e-6;

    // Vector of total disp. loadings
    std::vector<double> DLs = make_Fvec(disp0, delta_disp, steps_II);



    // (***************)
    // axial loading (zz)
    bool transverse_load = true;

    int steps_III = 9;     // steps of loading
    double tr_load0 = 1.0 * 1.0e-8;
    double delta_tr = 1.0 * 1.0e-8;

    // Vector of total tr. loadings
    std::vector<double> TLs = make_Fvec(tr_load0, delta_tr, steps_III);

    //parameter Tr_load;
    //parameter Tr_tot_load;



    // Material Properties
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // alpha_L in x,y,z : coeff of thermal expantion
    parameter alpha_th_0;
    alpha_th_0|Copper    = 17.7e-6;
    alpha_th_0|Silver    = 17e-6;
    alpha_th_0|Hastelloy = 13.4e-6;
    alpha_th_0|Buffer    = 9.5e-6;
    alpha_th_0|YBCO      = 11e-6;


    // E : elastic modulus
    parameter E;
    E|Copper    = 80e9;
    E|Silver    = 70e9; // 69-74
    E|Hastelloy = 200e9;  // Ref. ??
    E|Buffer    = 170e9;  // Ref. ??
    E|YBCO      = 157e9;  // Ref. ??


    // nu : Poisson’s ratio
    parameter nu;
    nu|Copper    = 0.338;
    nu|Silver    = 0.367; // 0.365-0.390
    nu|Hastelloy = 0.307;
    nu|Buffer    = 0.226; // 0.365-0.390
    nu|YBCO      = 0.3;   // Ref. ??


    // shear modulus
    expression mu;
    mu = E /2.0 /(1.0+nu);

    // tangent modulus
    parameter Et;
    Et|Copper    = 4e9;
    Et|Silver    = 0.0;
    Et|Hastelloy = 8.5e9;
    Et|Buffer    = 0.0;
    Et|YBCO      = 0.0;


    // hardening modulus (or generalized plastic modulus)
    // H = 2/3 * Hp (Hp : plastic modulus)
    expression H;
    H = 2.0/3.0 * E*Et/(E-Et);


    // hardening model's identifier
    // isotropic hardening beta = 1
    // kinematic hardening beta = 0
    parameter beta;
    beta|body = 1.0;

    // intial yield value
    parameter Sy0;
    Sy0|Copper    = 275e6;  // Yield Strength: 62-69 MPa //https://www.efunda.com/materials/alloys/copper/properties.cfm
    Sy0|Silver    = 140e6;  //110-340
    Sy0|Hastelloy = 980e6;  //https://elginfasteners.com/resources/materials/material-specifications/hastelloy-c-276/
    Sy0|Buffer    = 40e6;   //110-340
    Sy0|YBCO      = 40e6;   //Konstantina_Konstantopoulou_1_PhDThesis
    //Sy0|Buffer    = 1100e6;
    //Sy0|YBCO      = 1100e6;



    // C : elastic modulus tensor
    expression C(6,6, {1.0-nu,     nu,     nu,                 0,                0,                0,
                           nu, 1.0-nu,     nu,                 0,                0,                0,
                           nu,     nu, 1.0-nu,                 0,                0,                0,
                            0,      0,      0,  0.5*(1.0-2.0*nu),                0,                0,
                            0,      0,      0,                 0, 0.5*(1.0-2.0*nu),                0,
                            0,      0,      0,                 0,                0, 0.5*(1.0-2.0*nu)});

    C = E/(1.0+nu)/(1.0-2.0*nu)*C;


    // vector of coeff of thermal expantion
    expression alpha_th(6,1, { alpha_th_0,
                               alpha_th_0,
                               alpha_th_0,
                                0.0,
                                0.0,
                                0.0} );


    // Fields & Boundary conditions
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // x,y & z fields
    field x("x"), y("y"), z("z");

    // incremental dislacement field (for single load elastic problem and equilibrium)
    field du("h1xyz");
    // shape functions order
    du.setorder(body, 1);
    // BC : z clamping
    du.setconstraint(st);


    // dummy field for plasticity trigger formulation
    field u_dummy("h1");
    // shape functions order
    u_dummy.setorder(body, 1);
    // BC : z clamping
    u_dummy.setconstraint(st);

    // for equilibrium problem (residual & Jacobian)
    field u("h1xyz");
    // shape functions order
    u.setorder(body, 1);
    // BC : z clamping
    u.setconstraint(st);

    // initial u_tot
    //u_tot.write(body, "./vtu_folder/u1000.vtu", 1);

    // Expressions & "h1d" Fields
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // "h1d" fields for:
    // 1) strain_tot , 2) stress_corr , 3) strain_p
    // 4) alpha , 5) sigma_y , 6) stress_dummy , 7) CTO
    std::vector<field> strain_tot(6);     // total strain tensor
    std::vector<field> stress_corr(6);    // total corrected stress
    std::vector<field> strain_p(6);       // total plastic strain
    std::vector<field> alpha(6);          // back stress
    std::vector<field> kappa(1);          // yield stress

    std::vector<field> stress_dummy(6);   // total stress tensor for equ. form.
    std::vector<field> CTO(36);           // sub. of consistent tangent operator


    // setting order for "h1d" field and initializing them
    for (int i = 0; i < 6; i++) {

        //
        strain_tot[i] = field("h1d");
        strain_tot[i].setorder(body, 2);
        strain_tot[i].setvalue(body, 0.0);
        //
        stress_corr[i] = field("h1d");
        stress_corr[i].setorder(body, 2);
        stress_corr[i].setvalue(body, 0.0);
        //
        strain_p[i] = field("h1d");
        strain_p[i].setorder(body, 2);
        strain_p[i].setvalue(body, 0.0);
        //
        alpha[i] = field("h1d");
        alpha[i].setorder(body, 2);
        alpha[i].setvalue(body, 0.0);
        //
        stress_dummy[i] = field("h1d");
        stress_dummy[i].setorder(body, 2);
        stress_dummy[i].setvalue(body, 0.0);
    }

    kappa[0] = field("h1d");
    kappa[0].setorder(body, 2);
    kappa[0].setvalue(body, sqrt(2.0/3.0)*Sy0);

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {

            CTO[i*6+j] = field("h1d");
            CTO[i*6+j].setorder(body, 2);
            CTO[i*6+j].setvalue(body, 0.0);
            //CTO[i*6+j].setvalue(body, entry(i,j, C));

        }
    }



    // concatenating "h1d" fields for plasticity & sum_strain_stress functions
    std::vector<field> all_fields_plastic, all_fields_elastic;
    all_fields_plastic = concatenate({strain_tot, stress_corr, strain_p, alpha, kappa, stress_dummy, CTO});
    all_fields_elastic = concatenate({strain_tot, stress_corr});

    // exp(s) representation for "h1d" fields

    // total strain as exp
    expression strain_tot_exp(6,1, { strain_tot[0], strain_tot[1], strain_tot[2],
                                     strain_tot[3], strain_tot[4], strain_tot[5] });

    // total corrected stress as exp
    expression stress_corr_exp(6,1, { stress_corr[0], stress_corr[1], stress_corr[2],
                                      stress_corr[3], stress_corr[4], stress_corr[5] });

    // plastic strain as exp
    expression strain_p_exp(6,1, { strain_p[0], strain_p[1], strain_p[2],
                                   strain_p[3], strain_p[4], strain_p[5] });

    // back stress as exp
    expression alpha_exp(6,1, { alpha[0], alpha[1], alpha[2],
                                alpha[3], alpha[4], alpha[5] });

    // kappa as exp
    expression kappa_exp(1,1, {kappa[0]});

    // total stress tensor as exp
    expression stress_dummy_exp(6,1, { stress_dummy[0], stress_dummy[1], stress_dummy[2],
                                      stress_dummy[3], stress_dummy[4], stress_dummy[5] });
    // sub. of CTO as exp
    expression CTO_exp(6,6, { CTO[0],  CTO[1],  CTO[2],  CTO[3],  CTO[4],  CTO[5],
                              CTO[6],  CTO[7],  CTO[8],  CTO[9],  CTO[10], CTO[11],
                              CTO[12], CTO[13], CTO[14], CTO[15], CTO[16], CTO[17],
                              CTO[18], CTO[19], CTO[20], CTO[21], CTO[22], CTO[23],
                              CTO[24], CTO[25], CTO[26], CTO[27], CTO[28], CTO[29],
                              CTO[30], CTO[31], CTO[32], CTO[33], CTO[34], CTO[35] });



    // thermal normal stress to be subtracted
    parameter delta_sig_th(6,1);
    delta_sig_th|body = 1.0 * E / (1.0-2.0*nu) * alpha_th * delta_T;

    // single step elastic strain tensor
    expression strain_e = strain(du);

    // single step elastic stress tensor
    parameter stress_e(6,1);
    stress_e|body = C*strain(du) - delta_sig_th;

    // flag for updating p.box
    parameter pbox_update_params;
    pbox_update_params|body = 0;

    // the exp used in plasticity_trigger which ends in plasticity function
    expression exp_for_plasticity(1,1, plasticity,
                    {strain_e, stress_e, mu, H, beta,
                     strain_tot_exp, stress_corr_exp,
                     strain_p_exp, alpha_exp, kappa_exp, CTO_exp,
                     x, y, z, pbox_update_params},

                     all_fields_plastic);


    // the exp used in linear_loading_sum_up which ends in sum_strain_stress function
    expression exp_for_sum_up(1,1, sum_strain_stress,
                    {strain_e, stress_e,
                     strain_tot_exp, stress_corr_exp},

                     all_fields_elastic);

    // Formulation(s)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // cool-down
    // # # # # # # # # # # # # # #
    // # # # # # # # # # # # # # #
    // elastic weak formulation for single step cooling-down (from T0 to T)
    formulation elasticity_cooldown;
    elasticity_cooldown += integral(body, - C * strain(dof(du)) * strain(tf(du)));
    elasticity_cooldown += integral(body,  delta_sig_th * strain(tf(du)));

    // equilibrium formulation
    // *** By computing B (RHS), r(u_i) will be computed: r(u_i) ~ 0
    // *** A (stiffnessmatrix) can be used as the Jacobian for the newton method!
    formulation equilibrium_cooldown;
    // LHS: stiffnessmatrix as Jacobian
    equilibrium_cooldown += integral(body, - (C - CTO_exp) * strain(dof(u)) * strain(tf(u))); // A_as_J
    // RHS: residual
    equilibrium_cooldown += integral(body, - stress_dummy_exp * strain(tf(u)));                // B_of_eq

    // axial loading
    // # # # # # # # # # # # # # #
    // # # # # # # # # # # # # # #
    // single step elasticity problem for axial loading
    formulation elasticity_axial;
    elasticity_axial += integral(body, - C * strain(dof(du)) * strain(tf(du)) );           // A_dF

    // equilibrium formulation
    // *** By computing B (RHS), r(u_i) will be computed: r(u_i) ~ 0
    // *** A (stiffnessmatrix) can be used as the Jacobian for the newton method!
    formulation equilibrium;
    // LHS: stiffnessmatrix as Jacobian
    equilibrium += integral(body, - (C - CTO_exp) * strain(dof(u)) * strain(tf(u)) );   // A_as_J
    // RHS: residual
    equilibrium += integral(body,  - stress_dummy_exp * strain(tf(u)) );

    // transverse loading
    // # # # # # # # # # # # # # #
    // # # # # # # # # # # # # # #
    // single step elasticity problem for transverse loading
    formulation elasticity_transverse;
    elasticity_transverse += integral(body, - C * strain(dof(du)) * strain(tf(du)) );
    //elasticity_transverse += integral(mid_top,  - array1x3(0.0, Tr_load, 0.0)*tf(du));

    formulation equilibrium_transverse;
    // LHS: stiffnessmatrix as Jacobian
    equilibrium_transverse += integral(body, - (C - CTO_exp) * strain(dof(u)) * strain(tf(u)) );   // A_as_J
    // RHS: residual
    equilibrium_transverse += integral(body, - stress_dummy_exp * strain(tf(u)) );
    //equilibrium_transverse += integral(mid_top,  - array1x3(0.0, Tr_tot_load, 0.0)*tf(u));



    // dummy formulations for CF trigging
    // # # # # # # # # # # # # # #
    // # # # # # # # # # # # # # #
    // *** By computing B (RHS), the plasticity CF will be triggered
    formulation plasticity_trigger;
    plasticity_trigger += integral(body, exp_for_plasticity * tf(u_dummy));

    // *** By computing B (RHS), the sum_strain_stress CF will be triggered
    formulation linear_loading_sum_up;
    linear_loading_sum_up += integral(body, exp_for_sum_up * tf(u_dummy));


    // A(s) & B(s) & sol(s)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    // for "plasticity_trigger" & "linear_loading_sum_up"
    vec dummy_B;


    // for "elasticity_cooldown"
    mat A_cd;
    vec B_cd;

    vec sol_cd;
    vec sol_cd_final;

    // for "elasticity_axial"
    mat A_ax;
    vec B_ax;
    vec B_ax0;

    vec sol_ax0;
    vec sol_ax;

    // for "elasticity_transversal"
    mat A_tr;
    vec B_tr;
    vec B_tr0;

    vec sol_tr0;
    vec sol_tr;

    // for "equilibrium"
    mat A_as_J;
    vec B_of_eq;

    // for "equilibrium"
    vec source;


    // Loading (I)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (th_load == true) {

        //
        std::cout << "\n***Cool-down Loading***" << '\n';

        // delta_T for all incremental cooldown steps
        delta_T |body = -dT;

        // A & b >> sol.
        elasticity_cooldown.generatestiffnessmatrix();
        A_cd = elasticity_cooldown.A();

        elasticity_cooldown.generaterhs();
        B_cd = elasticity_cooldown.b();

        sol_cd = solve(A_cd, B_cd);


        // delta_T for the last step of cooldown
        delta_T |body = dT_final;

        // A & b >> sol.
        elasticity_cooldown.generatestiffnessmatrix();
        A_cd = elasticity_cooldown.A();

        elasticity_cooldown.generaterhs();
        B_cd = elasticity_cooldown.b();

        sol_cd_final = solve(A_cd, B_cd);



        switch(study_case) {

            case 0:
                //
                std::cout << "\nElastic Study:" << '\n';

                std::cout << "\ncooling down ..." << '\n';
                std::cout << "------------------" << '\n';

                // saving the zero state
                (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(100) +".vtu", 1);

                for (int i = 1 ; i<steps_I ; i++) {

                    // cooling from T_0 to T
                    std::cout << vecT[i-1]<< " K  >>>  " << vecT[i] << " K" << '\n';

                    if (i < steps_I - 1 ) {
                        //
                        delta_T |body = -dT;
                        du.setdata(body, sol_cd);

                        linear_loading_sum_up.generaterhs();
                        dummy_B = linear_loading_sum_up.b();

                    } else {

                        //
                        delta_T |body = dT_final;
                        du.setdata(body, sol_cd_final);

                        linear_loading_sum_up.generaterhs();
                        dummy_B = linear_loading_sum_up.b();
                    }


                    (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(100+i) +".vtu", 1);

                    extractdata(mymesh, body, avg_box, vol,
                                strain_tot_exp, stress_corr_exp,
                                thick_lays,
                                w, len);


                }
                std::cout << "Cool-down, Done!" << endl;
                break;





            case 1:

                //
                std::cout << "\nElastoplastic Study:" << '\n';

                // number of unknowns for nonlinear problem
                int num_of_u = B_cd.size();

                // petsc Vec(s) for disp.(s)
                // x0  & x : previous and current total disp : for loading loop   (i indexing)
                // y_0 & y : previous and current total disp : for nonlinear loop (j indexing)
                // du_acc  for nonlinear loop : (accumulated delta u)

                // zero value for initializing petsc Vec(s)
                PetscScalar ix0 = 0.0;

                // declaring x0 as a petsc Vec
                Vec x0;
                VecCreate(PETSC_COMM_WORLD,&x0);
                VecSetSizes(x0,PETSC_DECIDE,num_of_u);
                VecSetFromOptions(x0);

                // declaring other petsc Vec(s) using x0
                Vec x, y, y_0, du_acc;
                VecDuplicate(x0,&x);
                VecDuplicate(x0,&y);
                VecDuplicate(x0,&y_0);
                VecDuplicate(x0,&du_acc);


                // initializing both x & x0 with zero values
                VecSet(x0, ix0);
                VecSet(x, ix0);

                double loop_check;      // nonlinear loop exit check value
                bool NewRaph_flag;         // nonlinear loop exit flag
                int itr_num;            // iteration nums. of nonlinear loop
                PetscScalar res_norm_0; // residual norm of previous iteration in NL loop


                // cooling-down
                std::cout << "\ncooling down ..." << '\n';
                std::cout << "------------------" << '\n';

                // saving the zero state
                (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(100) +".vtu", 1);

                for (int i = 1 ; i<steps_I ; i++) {

                    // cooling from T_0 to T
                    std::cout << vecT[i-1]<< " K  >>>  " << vecT[i] << " K" << '\n';

                    // adding a linear (elastic) sol. to total disp.
                    if (i < steps_I - 1 ) {
                        // + initial linear sol. to x
                        VecAXPY(x, 1.0, sol_cd.getpetsc());
                        //
                        delta_T |body = -dT;
                    } else {
                        // + incremental linear sol. to x
                        VecAXPY(x, 1.0, sol_cd_final.getpetsc());
                        //
                        delta_T |body = dT_final;
                    }

                    // re-setting y_0 to x0 (for the NL_loop)
                    VecCopy(x0,y_0);
                    // re-setting du_acc to zero (for the NL_loop)
                    VecSet(du_acc, ix0);
                    // re-setting y to x (for the NL_loop)
                    VecCopy(x,y);

                    // nonlinear solving (as in SNESSolve(snes,NULL,x);)
                    itr_num = 0;
                    NewRaph_flag = true;

                    // re-setting the previous iteration res. norm to a large value
                    res_norm_0 = 1e10;

                    // NL loop
                    while (NewRaph_flag == true) {

                        // declaring local du_acc for Vec2vec (SL&petsc) problem!
                        Vec d_acc_loc;
                        VecDuplicate(x0,&d_acc_loc);
                        VecSet(d_acc_loc, ix0);

                        // making local du_acc with y, y_0 & real du_acc
                        VecAXPY(d_acc_loc, 1.0, y);
                        VecAXPY(d_acc_loc, -1.0, y_0);
                        VecAXPY(du_acc, 1, d_acc_loc);
                        VecCopy(du_acc, d_acc_loc);

                        // transforming local du_acc to SL vec (Vec2vec)
                        std::shared_ptr<rawvec> rv(new rawvec(elasticity_cooldown.getdofmanager(), d_acc_loc));
                        vec du_acc_vec(rv);

                        // setting du with accumulated delta u
                        du.setdata(body, du_acc_vec);

                        // plasticity trigger!
                        plasticity_trigger.generaterhs();
                        dummy_B = plasticity_trigger.b();
                        //abort();
                        // >>> new dummy stress made!
                        // >>> new CTO made!

                        // computing residual
                        equilibrium_cooldown.generaterhs();
                        B_of_eq = equilibrium_cooldown.b();

                        // computing the Jacobian with CTO made in plasticity
                        equilibrium_cooldown.generatestiffnessmatrix();
                        A_as_J = equilibrium_cooldown.A();

                        // A_as_J * B_of_eq :: J^-1 * r(u_j)
                        vec dy_vec;
                        dy_vec = solve(A_as_J, B_of_eq);

                        // updating y_0
                        VecCopy(y, y_0);
                        // new y
                        VecAXPY(y, 1.0, dy_vec.getpetsc());


                        // computing error
                        PetscScalar y0_norm, res_norm, J_res_norm;
                        VecNorm(y_0, NORM_2, &y0_norm);
                        VecNorm(B_of_eq.getpetsc(), NORM_2, &res_norm);
                        VecNorm(dy_vec.getpetsc(), NORM_2, &J_res_norm);

                        loop_check = J_res_norm/y0_norm;

                        printf("%*c Err_rel= %.2e   Err_res= %.2e\n", 25, ' ', loop_check, res_norm);
                        //printf("%.2e\n", res_norm);


                        // exit cases:

                        // 1) the current error is larger than previous iteration
                        // *********************
                        // *********************
                        if (res_norm_0 <= res_norm) {
                            //
                            NewRaph_flag = false;

                            // changing the disp. to the previous iteration
                            VecCopy(y_0, y);

                            std::cout << "\n!!! used the previous iteration ...\n" << '\n';
                        } else {
                            // saving the res. error of current iteration
                            res_norm_0 = res_norm;
                        }

                        // 2) the res. error reaches the threshold value
                        // *********************
                        // *********************
                        if (res_norm < 1.0e-7) {
                            //
                            NewRaph_flag = false;

                            // adding the last delta_u to du_acc to compute the stress by p.func.
                            if (itr_num > 0) {
                                //
                                Vec d_loc_;
                                VecDuplicate(x0,&d_loc_);
                                VecSet(d_loc_, ix0);

                                VecAXPY(d_loc_, 1.0, y);
                                VecAXPY(d_loc_, -1.0, y_0);
                                VecAXPY(du_acc, 1, d_loc_);
                                VecCopy(du_acc, d_loc_);

                                // transforming local du_acc to SL vec (Vec2vec)
                                std::shared_ptr<rawvec> rv_(new rawvec(elasticity_cooldown.getdofmanager(), d_loc_));
                                vec du_final(rv_);

                                // setting du with accumulated delta u
                                du.setdata(body, du_final);
                            }
                        }

                        // 3) the iterations reach the a max number!
                        // *********************
                        // *********************
                        if (itr_num > 200) {
                            //
                            NewRaph_flag = false;

                            std::cout << "\n\n    --- too many itrs! "<< '\n';
                        }

                        itr_num++;
                    }

                    // ********* one loading step done! *********
                    // ------------------------------------------
                    // ------------------------------------------

                    printf("%*c num_itr= %i\n", 4, ' ', itr_num-1);

                    // updating x
                    VecCopy(y, x);
                    // updating x0
                    VecCopy(x,x0);

                    // updating plasticity params.
                    pbox_update_params|body = 1;
                    plasticity_trigger.generaterhs();
                    dummy_B = plasticity_trigger.b();
                    pbox_update_params|body = 0;

                    // re-setting CTO to zero for next step of loading
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < 6; j++) {
                            CTO[i*6+j].setvalue(body, 0.0);
                        }
                    }


                    // writing u
                    //u.setdata(body, sol);

                    (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(100+i) +".vtu", 1);

                    extractdata(mymesh, body, avg_box, vol,
                                strain_tot_exp, stress_corr_exp,
                                thick_lays,
                                w, len);

                }

                std::cout << "Cool-down, Done!" << endl;
                break;

        }


    }


    // Loading (II)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (axial_disp == true) {

        //
        std::cout << "\n***Axial Loading***" << '\n';

        // re-defining stress for axial loading case
        stress_e|body = C*strain(du);

        // initial disp. loading
        du.compz().setconstraint(en, disp0);

        // A & b >> sol.
        elasticity_axial.generatestiffnessmatrix();
        A_ax = elasticity_axial.A();

        elasticity_axial.generaterhs();
        B_ax0 = elasticity_axial.b();

        sol_ax0 = solve(A_ax, B_ax0);


        // single step load
        du.compz().setconstraint(en, delta_disp);

        // A & b >> sol.
        elasticity_axial.generatestiffnessmatrix();
        A_ax = elasticity_axial.A();

        elasticity_axial.generaterhs();
        B_ax = elasticity_axial.b();

        sol_ax = solve(A_ax, B_ax);

        switch(study_case) {

            case 0:
                //
                std::cout << "\nElastic Study:" << '\n';

                std::cout << "\nLoading ..." << '\n';
                std::cout << "------------------" << '\n';

                // saving the zero state
                (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(100) +".vtu", 1);


                for (int i = 0; i<(steps_II+1); i++) {

                    // increasing the total load
                    vecDispLoad.push_back(DLs[i]);
                    printf("\nstep=%i   Disp. Loading : %.2e m\n",i, DLs[i]);

                    // adding a linear (elastic) sol. to total disp.
                    if (i == 0 ) {

                        du.setdata(body, sol_ax0);

                        linear_loading_sum_up.generaterhs();
                        dummy_B = linear_loading_sum_up.b();

                    } else {

                        du.setdata(body, sol_ax);

                        linear_loading_sum_up.generaterhs();
                        dummy_B = linear_loading_sum_up.b();

                    }


                    (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(101+i) +".vtu", 1);

                    extractdata(mymesh, body, avg_box, vol,
                                strain_tot_exp, stress_corr_exp,
                                thick_lays,
                                w, len);

                }
                std::cout << "Axial loading, Done!" << endl;
                break;

            case 1:

                //
                std::cout << "\nElastoplastic Study:" << '\n';

                // number of unknowns for nonlinear problem
                int num_of_u = B_ax.size();


                // petsc Vec(s) for disp.(s)
                // x0  & x : previous and current total disp : for loading loop   (i indexing)
                // y_0 & y : previous and current total disp : for nonlinear loop (j indexing)
                // du_acc  for nonlinear loop : (accumulated delta u)

                // zero value for initializing petsc Vec(s)
                PetscScalar ix0 = 0.0;

                // declaring x0 as a petsc Vec
                Vec x0;
                VecCreate(PETSC_COMM_WORLD,&x0);
                VecSetSizes(x0,PETSC_DECIDE,num_of_u);
                VecSetFromOptions(x0);

                // declaring other petsc Vec(s) using x0
                Vec x, y, y_0, du_acc;
                VecDuplicate(x0,&x);
                VecDuplicate(x0,&y);
                VecDuplicate(x0,&y_0);
                VecDuplicate(x0,&du_acc);


                // initializing both x & x0 with zero values
                VecSet(x0, ix0);
                VecSet(x, ix0);

                double loop_check;      // nonlinear loop exit check value
                bool NewRaph_flag;         // nonlinear loop exit flag
                int itr_num;            // iteration nums. of nonlinear loop
                PetscScalar res_norm_0; // residual norm of previous iteration in NL loop

                // loading
                std::cout << "\n\nLoading ..." << '\n';
                std::cout << "------------------" << '\n';

                // saving the zero state
                (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(100) +".vtu", 1);

                for (int i = 0; i<(steps_II+1); i++) {

                    // increasing the total load
                    vecDispLoad.push_back(DLs[i]);
                    printf("\nstep=%i   Disp. Loading : %.2e m\n",i, DLs[i]);

                    // set the disp load for equilibrium (residual) formualation and Jacobian formulation
                    u.compz().setconstraint(en, DLs[i]);

                    // adding a linear (elastic) sol. to total disp.
                    if (i == 0 ) {
                        // + initial linear sol. to x
                        VecAXPY(x, 1.0, sol_ax0.getpetsc());

                        // initializing the source vector by RHS of intial disp. load
                        source = B_ax;

                    } else {
                        // + incremental linear sol. to x
                        VecAXPY(x, 1.0, sol_ax.getpetsc());

                        // + source vector by RHS of incremental disp. load
                        source = source + B_ax;
                    }

                    // re-setting y_0 to x0 (for the NL_loop)
                    VecCopy(x0,y_0);
                    // re-setting du_acc to zero (for the NL_loop)
                    VecSet(du_acc, ix0);
                    // re-setting y to x (for the NL_loop)
                    VecCopy(x,y);


                    // nonlinear solving (as in SNESSolve(snes,NULL,x);)
                    itr_num = 0;
                    NewRaph_flag = true;

                    // re-setting the previous iteration res. norm to a large value
                    res_norm_0 = 1e6;

                    // NL loop
                    while (NewRaph_flag == true) {

                        // declaring local du_acc for Vec2vec (SL&petsc) problem!
                        Vec d_acc_loc;
                        VecDuplicate(x0,&d_acc_loc);
                        VecSet(d_acc_loc, ix0);

                        // making local du_acc with y, y_0 & real du_acc
                        VecAXPY(d_acc_loc, 1.0, y);
                        VecAXPY(d_acc_loc, -1.0, y_0);
                        VecAXPY(du_acc, 1, d_acc_loc);
                        VecCopy(du_acc, d_acc_loc);

                        // transforming local du_acc to SL vec (Vec2vec)
                        std::shared_ptr<rawvec> rv(new rawvec(elasticity_axial.getdofmanager(), d_acc_loc));
                        vec du_acc_vec(rv);

                        // setting du with accumulated delta u
                        du.setdata(body, du_acc_vec);

                        // plasticity trigger!
                        plasticity_trigger.generaterhs();
                        dummy_B = plasticity_trigger.b();
                        // >>> new dummy stress made!
                        // >>> new CTO made!

                        // computing residual
                        equilibrium.generaterhs();
                        B_of_eq = equilibrium.b();

                        // adding the source term to the residual
                        B_of_eq = B_of_eq - source;

                        // computing the Jacobian with CTO made in plasticity
                        equilibrium.generatestiffnessmatrix();
                        A_as_J = equilibrium.A();

                        // A_as_J * B_of_eq :: J^-1 * r(u_j)
                        vec dy_vec;
                        dy_vec = solve(A_as_J, B_of_eq);

                        // updating y_0
                        VecCopy(y, y_0);
                        // new y
                        VecAXPY(y, 1.0, dy_vec.getpetsc());


                        // computing error
                        PetscScalar y0_norm, res_norm, J_res_norm;
                        VecNorm(y_0, NORM_2, &y0_norm);
                        VecNorm(B_of_eq.getpetsc(), NORM_2, &res_norm);
                        VecNorm(dy_vec.getpetsc(), NORM_2, &J_res_norm);

                        loop_check = J_res_norm/y0_norm;

                        printf("%*c Err_rel= %.2e   Err_res= %.2e\n", 25, ' ', loop_check, res_norm);
                        //printf("%.2e\n", res_norm);


                        // exit cases:

                        // 1) the current error is larger than previous iteration
                        // *********************
                        // *********************
                        if (res_norm_0 <= res_norm) {
                            //
                            NewRaph_flag = false;

                            // changing the disp. to the previous iteration
                            VecCopy(y_0, y);

                            std::cout << "\n!!! used the previous iteration ...\n" << '\n';
                        } else {
                            // saving the res. error of current iteration
                            res_norm_0 = res_norm;
                        }

                        // 2) the res. error reaches the threshold value
                        // *********************
                        // *********************
                        if (res_norm < 1.0e-7) {
                            //
                            NewRaph_flag = false;

                            // adding the last delta_u to du_acc to compute the stress by p.func.
                            if (itr_num > 0) {
                                //
                                Vec d_loc_;
                                VecDuplicate(x0,&d_loc_);
                                VecSet(d_loc_, ix0);

                                VecAXPY(d_loc_, 1.0, y);
                                VecAXPY(d_loc_, -1.0, y_0);
                                VecAXPY(du_acc, 1, d_loc_);
                                VecCopy(du_acc, d_loc_);

                                // transforming local du_acc to SL vec (Vec2vec)
                                std::shared_ptr<rawvec> rv_(new rawvec(elasticity_axial.getdofmanager(), d_loc_));
                                vec du_final(rv_);

                                // setting du with accumulated delta u
                                du.setdata(body, du_final);
                            }
                        }

                        // 3) the iterations reach the a max number!
                        // *********************
                        // *********************
                        if (itr_num > 200) {
                            //
                            NewRaph_flag = false;

                            std::cout << "\n\n    --- too many itrs! "<< '\n';
                        }

                        itr_num++;
                    }

                    // ******** one  loading step done! *********
                    // ------------------------------------------
                    // ------------------------------------------

                    printf("%*c num_itr= %i\n", 4, ' ', itr_num-1);

                    // updating x
                    VecCopy(y, x);
                    // updating x0
                    VecCopy(x,x0);

                    // updating plasticity params.
                    pbox_update_params|body = 1;
                    plasticity_trigger.generaterhs();
                    dummy_B = plasticity_trigger.b();
                    pbox_update_params|body = 0;

                    // re-setting CTO to zero for next step of loading
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < 6; j++) {
                            CTO[i*6+j].setvalue(body, 0.0);
                        }
                    }

                    (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(101+i) +".vtu", 1);

                    extractdata(mymesh, body, avg_box, vol,
                                strain_tot_exp, stress_corr_exp,
                                thick_lays,
                                w, len);


                }
                std::cout << "Axial loading, Done!" << endl;
                break;

        }

    }


    // Loading (III)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (transverse_load == true) {

        //
        std::cout << "\n***Transversal Loading***" << '\n';

        // B.C.
        // * # * # * # * # * # * # * # * #
        // * # * # * # * # * # * # * # * #
        du.setconstraint(en);

        double y_bot_limit = 0;
        expression bot_condexpr_du = ifpositive(-compy(du), 1, -1);
        du.compy().setconditionalconstraint(bot, bot_condexpr_du, y_bot_limit-1e-9);

        u.setconstraint(en);

        expression bot_condexpr_u = ifpositive(-compy(u), 1, -1);
        u.compy().setconditionalconstraint(bot, bot_condexpr_u, y_bot_limit-1e-9);


        // re-defining stress
        // * # * # * # * # * # * # * # * #
        // * # * # * # * # * # * # * # * #
        stress_e|body = C*strain(du);


        // intial load
        // * # * # * # * # * # * # * # * #
        // * # * # * # * # * # * # * # * #
        du.compy().setconstraint(mid_top, -tr_load0);

        // A & b >> sol.
        elasticity_transverse.generatestiffnessmatrix();
        A_tr = elasticity_transverse.A();

        elasticity_transverse.generaterhs();
        B_tr0 = elasticity_transverse.b();

        sol_tr0 = solve(A_tr, B_tr0);


        // single step load
        // * # * # * # * # * # * # * # * #
        // * # * # * # * # * # * # * # * #
        du.compy().setconstraint(mid_top, -delta_tr);


        // A & b >> sol.
        elasticity_transverse.generatestiffnessmatrix();
        A_tr = elasticity_transverse.A();

        elasticity_transverse.generaterhs();
        B_tr = elasticity_transverse.b();

        sol_tr = solve(A_tr, B_tr);


        // saving the zero state
        //(entry(0,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/sxx"+ std::to_string(100) +".vtu", 1);
        //(entry(1,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/syy"+ std::to_string(100) +".vtu", 1);
        //(entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(100) +".vtu", 1);
        //du.write(body,"./vtu_folder/u"+ std::to_string(100) +".vtu", 1);

        switch(study_case) {

            case 0:
                //
                std::cout << "\nElastic Study:" << '\n';

                std::cout << "\nLoading ..." << '\n';
                std::cout << "------------------" << '\n';

                for (int i = 0 ; i<(steps_III+1) ; i++) {

                    // increasing the total load
                    vecTrForce.push_back(TLs[i]);
                    printf("\nstep=%i   Tr. disp. Loading : %1.3e m\n",i, TLs[i]);

                    du.compy().setconstraint(mid_top, -TLs[i]);

                    // A & b >> sol.
                    elasticity_transverse.generatestiffnessmatrix();
                    A_tr = elasticity_transverse.A();

                    elasticity_transverse.generaterhs();
                    B_tr0 = elasticity_transverse.b();

                    sol_tr0 = solve(A_tr, B_tr0);

                    du.setdata(body, sol_tr0);

                    du.write(body,"./vtu_folder/u"+ std::to_string(101+i) +".vtu", 1);
                    //(entry(0,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/sxx"+ std::to_string(101+i) +".vtu", 1);
                    //(entry(1,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/syy"+ std::to_string(101+i) +".vtu", 1);
                    //(entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(101+i) +".vtu", 1);

                    /*
                    extractdata(mymesh, body, mid_box, vol_mid_box,
                                strain_tot_exp, stress_corr_exp,
                                thick_lays,
                                w, len);*/


                    //

                }
                std::cout << "Transverse loading, Done!" << endl;
                break;

            case 1:

                //
                std::cout << "\nElastoplastic Study:" << '\n';

                // number of unknowns for nonlinear problem
                int num_of_u = B_tr.size();


                // petsc Vec(s) for disp.(s)
                // x0  & x : previous and current total disp : for loading loop   (i indexing)
                // y_0 & y : previous and current total disp : for nonlinear loop (j indexing)
                // du_acc  for nonlinear loop : (accumulated delta u)

                // zero value for initializing petsc Vec(s)
                PetscScalar ix0 = 0.0;

                // declaring x0 as a petsc Vec
                Vec x0;
                VecCreate(PETSC_COMM_WORLD,&x0);
                VecSetSizes(x0,PETSC_DECIDE,num_of_u);
                VecSetFromOptions(x0);

                // declaring other petsc Vec(s) using x0
                Vec x, y, y_0, du_acc;
                VecDuplicate(x0,&x);
                VecDuplicate(x0,&y);
                VecDuplicate(x0,&y_0);
                VecDuplicate(x0,&du_acc);


                // initializing both x & x0 with zero values
                VecSet(x0, ix0);
                VecSet(x, ix0);

                double loop_check;        // nonlinear loop exit check value
                bool NewRaph_flag;        // nonlinear loop exit flag
                int itr_num;              // iteration nums. of nonlinear loop
                PetscScalar res_norm_0;   // residual norm of previous iteration in NL loop
                bool initial_flag = true; // flag for: initial load or delta load
                bool save_flag;           // this flag allows to reset the initial load
                bool second_restart = true;

                // loading
                std::cout << "\n\nLoading ..." << '\n';
                std::cout << "------------------" << '\n';

                for (int i = 0; i<(steps_III+1); i++) {

                    // increasing the total load
                    vecTrForce.push_back(TLs[i]);

                    // total trans. load for eq. formulation
                    u.compy().setconstraint(mid_top, -TLs[i]);


                    // assigning true value to continue
                    save_flag = true;

                    // no indenting for save_flag loop
                    while (save_flag) {
                    //
                    printf("\nstep=%i   Tr. disp. Loading : %1.3e m\n",i, TLs[i]);

                    // adding a linear (elastic) sol. to total sol.
                    if (initial_flag) {
                        // + initial linear sol. to x
                        VecAXPY(x, 1.0, sol_tr0.getpetsc());

                        // initializing the source vector by RHS of intial disp. load
                        source = B_tr0;

                        // changing the flag for next step to use delta load
                        initial_flag = false;

                    } else {
                        // + incremental linear sol. to x
                        VecAXPY(x, 1.0, sol_tr.getpetsc());

                        // + source vector by RHS of incremental disp. load
                        source = source + B_tr;
                    }

                    // re-setting y_0 to x0 (for the NL_loop)
                    VecCopy(x0,y_0);
                    // re-setting du_acc to zero (for the NL_loop)
                    VecSet(du_acc, ix0);
                    // re-setting y to x (for the NL_loop)
                    VecCopy(x,y);


                    // nonlinear solving (as in SNESSolve(snes,NULL,x);)
                    itr_num = 0;
                    NewRaph_flag = true;

                    // re-setting the previous iteration res. norm to a large value
                    res_norm_0 = 1e10;

                    // NL loop
                    while (NewRaph_flag == true) {

                        // declaring local du_acc for Vec2vec (SL&petsc) problem!
                        Vec d_acc_loc;
                        VecDuplicate(x0,&d_acc_loc);
                        VecSet(d_acc_loc, ix0);

                        // making local du_acc with y, y_0 & real du_acc
                        VecAXPY(d_acc_loc, 1.0, y);
                        VecAXPY(d_acc_loc, -1.0, y_0);
                        VecAXPY(du_acc, 1, d_acc_loc);
                        VecCopy(du_acc, d_acc_loc);

                        // transforming local du_acc to SL vec (Vec2vec)
                        std::shared_ptr<rawvec> rv(new rawvec(elasticity_transverse.getdofmanager(), d_acc_loc));
                        vec du_acc_vec(rv);

                        // setting du with accumulated delta u
                        du.setdata(body, du_acc_vec);


                        // plasticity trigger!
                        plasticity_trigger.generaterhs();
                        dummy_B = plasticity_trigger.b();
                        // >>> new dummy stress made!
                        // >>> new CTO made!

                        // computing residual
                        equilibrium_transverse.generaterhs();
                        B_of_eq = equilibrium_transverse.b();

                        // adding the source term to the residual
                        B_of_eq = B_of_eq - source;

                        // computing the Jacobian with CTO made in plasticity
                        equilibrium_transverse.generatestiffnessmatrix();
                        A_as_J = equilibrium_transverse.A();

                        // A_as_J * B_of_eq :: J^-1 * r(u_j)
                        vec dy_vec;
                        dy_vec = solve(A_as_J, B_of_eq);

                        // updating y_0
                        VecCopy(y, y_0);
                        // new y
                        VecAXPY(y, 1.0, dy_vec.getpetsc());


                        // computing error
                        PetscScalar y0_norm, res_norm, J_res_norm;
                        VecNorm(y_0, NORM_2, &y0_norm);
                        VecNorm(B_of_eq.getpetsc(), NORM_2, &res_norm);
                        VecNorm(dy_vec.getpetsc(), NORM_2, &J_res_norm);

                        loop_check = J_res_norm/y0_norm;

                        printf("%*c Err_rel= %.2e   Err_res= %.2e\n", 25, ' ', loop_check, res_norm);
                        //printf("%.2e\n", res_norm);


                        // exit cases:

                        // 1) the current error is larger than previous iteration
                        // *********************
                        // *********************
                        if (res_norm_0 <= res_norm) {

                            //
                            std::cout << "\n!!! stopped for not converging ...\n" << '\n';

                            if (second_restart) {
                                // changing flags for exit both loops
                                NewRaph_flag = false;
                                save_flag = false;

                                //
                                second_restart = false;

                                std::cout << "\nThe initial load is already restarted! It's skipped!" << '\n';
                            } else {
                                //
                                // changing flags
                                NewRaph_flag = false; // to exit NewRaph loop (go back to the start of save_flag)
                                initial_flag = true;  // to resest the initial load

                                // resetting all vecs
                                // some of them is not needed! >>> just to be sure!
                                VecSet(x0, ix0);
                                VecSet(x, ix0);
                                VecSet(y_0, ix0);
                                VecSet(y, ix0);
                                VecSet(du_acc, ix0);

                                // resetting du
                                Vec zero_Vec;
                                VecDuplicate(x0,&zero_Vec);
                                VecSet(zero_Vec, ix0);
                                std::shared_ptr<rawvec> rv_0(new rawvec(elasticity_transverse.getdofmanager(), zero_Vec));
                                vec zero_vec_(rv_0);
                                du.setdata(body, zero_vec_);

                                // solving a new initial load
                                du.compy().setconstraint(mid_top, -TLs[i]);

                                // A & b >> sol.
                                elasticity_transverse.generatestiffnessmatrix();
                                A_tr = elasticity_transverse.A();

                                elasticity_transverse.generaterhs();
                                B_tr0 = elasticity_transverse.b();

                                sol_tr0 = solve(A_tr, B_tr0);

                                // resetting all field values
                                for (int jj = 0; jj < 6; jj++) {
                                    //
                                    strain_tot[jj].setvalue(body, 0.0);
                                    stress_corr[jj].setvalue(body, 0.0);
                                    strain_p[jj].setvalue(body, 0.0);
                                    alpha[jj].setvalue(body, 0.0);
                                    stress_dummy[jj].setvalue(body, 0.0);
                                    //
                                }

                                kappa[0].setvalue(body, sqrt(2.0/3.0)*Sy0);

                                for (int ii = 0; ii < 6; ii++) {
                                    for (int jj = 0; jj < 6; jj++) {
                                        CTO[ii*6+jj].setvalue(body, 0.0);
                                    }
                                }

                                // using plasticity_trigger >>> just to be sure!
                                pbox_update_params|body = 1;
                                plasticity_trigger.generaterhs();
                                dummy_B = plasticity_trigger.b();
                                pbox_update_params|body = 0;

                                std::cout << "Same load as initial load!\n" << '\n';

                                //
                                second_restart = true;

                            }

                        } else {
                            // saving the res. error of current iteration
                            res_norm_0 = res_norm;

                        }

                        // 2) the res. error reaches the threshold value
                        // *********************
                        // *********************
                        if (res_norm < 1.0e-7) {
                            //
                            NewRaph_flag = false;
                            save_flag = false;

                            //
                            second_restart = false;

                            // adding the last delta_u to du_acc to compute the stress by p.func.
                            if (itr_num > 0) {
                                //
                                Vec d_loc_;
                                VecDuplicate(x0,&d_loc_);
                                VecSet(d_loc_, ix0);

                                VecAXPY(d_loc_, 1.0, y);
                                VecAXPY(d_loc_, -1.0, y_0);
                                VecAXPY(du_acc, 1, d_loc_);
                                VecCopy(du_acc, d_loc_);

                                // transforming local du_acc to SL vec (Vec2vec)
                                std::shared_ptr<rawvec> rv_(new rawvec(elasticity_transverse.getdofmanager(), d_loc_));
                                vec du_final(rv_);

                                // setting du with accumulated delta u
                                du.setdata(body, du_final);
                            }
                        }

                        // 3) the iterations reach the a max number!
                        // *********************
                        // *********************
                        if (itr_num > 200) {
                            //
                            NewRaph_flag = false;
                            save_flag = false;

                            std::cout << "\n\n    --- too many itrs! "<< '\n';
                        }

                        itr_num++;
                    }


                    } // save_flag : end

                    // ******** one  loading step done! *********
                    // ------------------------------------------
                    // ------------------------------------------

                    printf("%*c num_itr= %i\n", 4, ' ', itr_num-1);

                    // updating x
                    VecCopy(y, x);
                    // updating x0
                    VecCopy(x, x0);

                    // updating plasticity params.
                    pbox_update_params|body = 1;
                    plasticity_trigger.generaterhs();
                    dummy_B = plasticity_trigger.b();
                    pbox_update_params|body = 0;

                    // re-setting CTO to zero for next step of loading
                    for (int ii = 0; ii < 6; ii++) {
                        for (int jj = 0; jj < 6; jj++) {
                            CTO[ii*6+jj].setvalue(body, 0.0);
                        }
                    }

                    double val;
                    val = entry(2,0, strain_tot_exp).integrate(mid_box, 4);
                    //std::cout << "e_zz: " << val/vol_mid_box*100 <<  " %" << '\n';
                    //abort();


                    //expression my_von = vonmises(stress_corr_exp);
                    //double ybco_point = thCu + thAg + thHastelloy + thBuffer + thYBCO/2.0;
                    /*
                    val = entry(0,0, stress_corr_exp).interpolate(body, {w/2, ybco_point, len/2})[0];
                    std::cout << "s_xx: " << val/1e6 <<  " MPa" << '\n';
                    val = entry(1,0, stress_corr_exp).interpolate(body, {w/2, ybco_point, len/2})[0];
                    std::cout << "s_yy: " << val/1e6 <<  " MPa" << '\n';
                    val = entry(2,0, stress_corr_exp).interpolate(body, {w/2, ybco_point, len/2})[0];
                    std::cout << "s_zz: " << val/1e6 <<  " MPa" << '\n';

                    val = my_von.interpolate(body, {w/2, ybco_point, len/2})[0];
                    std::cout << "von: " << val/1e6 <<  " MPa" << '\n';

                    std::cout << " vol of mid box: " << vol_mid_box << '\n';

                    val = (entry(0,0, stress_corr_exp)/1e6).integrate(mid_box, 4);
                    std::cout << "s_xx: " << val/vol_mid_box <<  " MPa" << '\n';

                    val = (entry(1,0, stress_corr_exp)/1e6).integrate(mid_box, 4);
                    std::cout << "s_yy: " << val/vol_mid_box <<  " MPa" << '\n';

                    val = (entry(2,0, stress_corr_exp)/1e6).integrate(mid_box, 4);
                    std::cout << "s_zz: " << val/vol_mid_box <<  " MPa" << '\n';

                    lll = expression(1).integrate(cent_line, 4);

                    std::cout << " length_of_cent_line: " << lll << '\n';

                    val = entry(0,0, stress_corr_exp).integrate(cent_line, 4);
                    std::cout << "s_xx: " << val/lll <<  " MPa" << '\n';

                    val = entry(1,0, stress_corr_exp).integrate(cent_line, 4);
                    std::cout << "s_yy: " << val/lll <<  " MPa" << '\n';

                    val = entry(2,0, stress_corr_exp).integrate(cent_line, 4);
                    std::cout << "s_zz: " << val/lll <<  " MPa" << '\n';
                    */


                    (entry(0,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/sxx"+ std::to_string(101+i) +".vtu", 1);
                    (entry(1,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/syy"+ std::to_string(101+i) +".vtu", 1);
                    (entry(2,0, stress_corr_exp)/1e6).write(body,"./vtu_folder/szz"+ std::to_string(101+i) +".vtu", 1);



                    extractdata(mymesh, body, mid_box, vol_mid_box,
                                strain_tot_exp, stress_corr_exp,
                                thick_lays,
                                w, len);


                }

                std::cout << "Transverse loading, Done!" << endl;
                break;
        }
    }


    savedata(tags);

    myclock.print("\nTime elapsed: ");
    //abort();

}
// *****************************************************************************
// *****************************************************************************
mesh createmesh(double w, double len,
                vector<double> thick_lays,
                int nx, int nz, vector<int> num_nodes_lays)
{

    int Copper = 1, Silver = 2, Hastelloy = 3, Buffer = 4, YBCO = 5; // materials
    int body = 6;           // body : all tapes
    int st = 7, en = 8;     // tape's vertical surface at either sides
    int bot = 9, top = 10;  // bottom & top surfaces of the tape
    int mid_top = 11;       // top surface region for transverse force
    int avg_box = 12;       // vol. for avg. tape
    int cent_line = 13;     // center line for avg. behavior for transverse loading
    int mid_box = 14;       // the box subjected to the trans. load



    int ny;        // nx, ny, nz are the number of nodes in x, y, z direction in that shape
    double ya, yb; // the heights of each layer (from ya to yb)
    double xa, xb; // the x values for each layer (from xa to xb)

    // number of lay(s)
    int num_lays = thick_lays.size();

    // calculating thickAll
    double thickAll = 0.0;
    for (int i = 0; i<num_lays; i++) {
        thickAll = thick_lays[i] + thickAll;
    }

    // declaring vector of shapes for each quadrangle of the layers
    std::vector<shape> vecLs(num_lays);

    // generating the shapes based on coordinates & number of nodes
    xa = 0.0; xa = w;
    ya = 0.0; yb = 0;

    for (int i = 0 ; i<thick_lays.size() ; i++) {
        //
        ya = yb;
        yb = ya + thick_lays[i];
        // number of nodes
        ny = num_nodes_lays[i];
        //
        shape q("quadrangle", -1, {xa,ya,0, xb,ya,0, xb,yb,0, xa,yb,0}, {nx,ny,nx,ny});
        vecLs[i] = q.extrude(-1, len, nz);
    }

    // volumes as shape
    /*
    shape vCopper("union",    Copper,    {vecLs[0], vecLs[3]});
    shape vHastelloy("union", Hastelloy, {vecLs[1]});
    shape vYBCO("union",      YBCO,      {vecLs[2]});
    */

    shape vCopper("union",    Copper,    {vecLs[0], vecLs[6]});
    shape vSilver("union",    Silver,    {vecLs[1], vecLs[5]});
    shape vHastelloy("union", Hastelloy, {vecLs[2]});
    shape vBuffer("union",    Buffer,    {vecLs[3]});
    shape vYBCO("union",      YBCO,      {vecLs[4]});

    // whole body Shape (the tape)
    shape allvols("union", body, vecLs);

    mesh mymesh; // mesh obj created

    double eps = 1e-7;
    // st , en
    mymesh.selectbox(st, body, 2, {0.0-eps, w+eps,   0.0-eps, thickAll+eps,   0.0-eps, 0.0+eps});
    mymesh.selectbox(en, body, 2, {0.0-eps, w+eps,   0.0-eps, thickAll+eps,   len-eps, len+eps});
    //mid_top
    mymesh.selectbox(mid_top, body, 2, {0.0-eps, w+eps,   thickAll-eps, thickAll+eps,   len/2.0-w/2.0, len/2.0+w/2.0});
    // bottom
    mymesh.selectbox(bot, body, 2, {0.0-eps, w+eps,   0.0-eps, 0.0+eps,  0.0-eps, len+eps});
    // top
    mymesh.selectbox(top, body, 2, {0.0-eps, w+eps,   thickAll-eps, thickAll+eps,  0.0-eps, len+eps});
    // average box
    mymesh.selectbox(avg_box, body, 3, {0.0-eps, w+eps,   0.0-eps, thickAll+eps,  0.0+0.15*len, 0.0+0.85*len});
    // center line
    mymesh.selectbox(cent_line, body, 1, {w/2.0-10*eps, w/2.0+10*eps,   0.0-eps, thickAll+eps,  len/2.0-10*eps, len/2.0+10*eps});
    // mid box
    mymesh.selectbox(mid_box, body, 3, {0.0-eps, w+eps,   0.0-eps, thickAll+eps,  len/2.0-4.0e-4-eps, len/2.0+4.0e-4+eps});
    //mymesh.selectbox(mid_box, body, 2, {0.0-eps, w+eps,   0.0-eps, thickAll+eps,  len/2.0-eps, len/2.0+eps});


    //
    //mymesh.load({allvols, vCopper, vHastelloy, vYBCO});
    mymesh.load({allvols, vCopper, vSilver, vHastelloy, vBuffer, vYBCO});
    mymesh.write("HTStape.msh");

    return mymesh;
}
// *****************************************************************************
// *****************************************************************************
std::vector<field> concatenate(std::vector<std::vector<field>> input)
{
    // concatenate fields in one field vector

    int totlen = 0;
    for (int i = 0; i < input.size(); i++) {
        totlen += input[i].size();
    }

    std::vector<field> output(totlen);

    int index = 0;
    for (int i = 0; i < input.size(); i++)
    {
        for (int j = 0; j < input[i].size(); j++)
        {
            output[index] = input[i][j];
            index++;
        }
    }

    return output;
}
// *****************************************************************************
// *****************************************************************************
std::vector<double> make_Tvec(double T_0, double T_f, double dT)
{
    // this function generates the temperature cooling-down vector
    // based on T_0 , T_f & dT

    // temperature vec to be returned
    std::vector<double> Ts;

    do {
        //
        if (T_0 <= T_f) {
            Ts.push_back(T_f);

            return Ts;
        }

        Ts.push_back(T_0);

        T_0 -= dT;
        //
    } while (1);

}
// *****************************************************************************
// *****************************************************************************
std::vector<double> make_Fvec(double d0, double dd, int steps)
{
    // used for : 1) axial disp. load 2) transverse load

    // d0 : intial
    // dd : increment
    // steps

    std::vector<double> Fs;
    Fs.push_back(d0);

    for (int i = 0; i < steps; i++) {
        //
        Fs.push_back(d0 + (i+1)*dd);
    }
    return Fs;
}
// *****************************************************************************
// *****************************************************************************
void prepvecs(int num_lays)
{
    // this function prepares these vectors (for storing pointers):
    // 1) strain : strain_xx, strain_yy, strain_zz, strain_yz, strain_xz, strain_xy
    // 2) stress : stress_xx, stress_yy, stress_zz, stress_yz, stress_xz, stress_xy
    // 3) plastic strain : ep_xx, ep_yy, ep_zz, ep_yz, ep_xz, ep_xy
    // 4) von-Mises stress
    // 5) yield stress

    for (int i=0; i<num_lays ; i++) {

        // initializing new pointers for each layer

        L_strain_xx = new vector<double>;
        L_strain_yy = new vector<double>;
        L_strain_zz = new vector<double>;
        L_strain_yz = new vector<double>;
        L_strain_xz = new vector<double>;
        L_strain_xy = new vector<double>;

        L_stress_xx = new vector<double>;
        L_stress_yy = new vector<double>;
        L_stress_zz = new vector<double>;
        L_stress_yz = new vector<double>;
        L_stress_xz = new vector<double>;
        L_stress_xy = new vector<double>;

        // pushing back the new poiners to the "main vector of pointers"

        strain_xx.push_back(L_strain_xx);
        strain_yy.push_back(L_strain_yy);
        strain_zz.push_back(L_strain_zz);
        strain_yz.push_back(L_strain_yz);
        strain_xz.push_back(L_strain_xz);
        strain_xy.push_back(L_strain_xy);

        stress_xx.push_back(L_stress_xx);
        stress_yy.push_back(L_stress_yy);
        stress_zz.push_back(L_stress_zz);
        stress_yz.push_back(L_stress_yz);
        stress_xz.push_back(L_stress_xz);
        stress_xy.push_back(L_stress_xy);

      }

      // inserting the zero values to all vectors

      // in all layers
      for (int i=0; i<num_lays ; i++) {

          strain_xx[i]->insert(strain_xx[i]->begin(), 0);
          strain_yy[i]->insert(strain_yy[i]->begin(), 0);
          strain_zz[i]->insert(strain_zz[i]->begin(), 0);
          strain_yz[i]->insert(strain_yz[i]->begin(), 0);
          strain_xz[i]->insert(strain_xz[i]->begin(), 0);
          strain_xy[i]->insert(strain_xy[i]->begin(), 0);

          stress_xx[i]->insert(stress_xx[i]->begin(), 0);
          stress_yy[i]->insert(stress_yy[i]->begin(), 0);
          stress_zz[i]->insert(stress_zz[i]->begin(), 0);
          stress_yz[i]->insert(stress_yz[i]->begin(), 0);
          stress_xz[i]->insert(stress_xz[i]->begin(), 0);
          stress_xy[i]->insert(stress_xy[i]->begin(), 0);

      }


      // inserting zero values for whole tape data
      tape_strain_xx.insert(tape_strain_xx.begin(), 0);
      tape_strain_yy.insert(tape_strain_yy.begin(), 0);
      tape_strain_zz.insert(tape_strain_zz.begin(), 0);
      tape_strain_yz.insert(tape_strain_yz.begin(), 0);
      tape_strain_xz.insert(tape_strain_xz.begin(), 0);
      tape_strain_xy.insert(tape_strain_xy.begin(), 0);
      //
      tape_stress_xx.insert(tape_stress_xx.begin(), 0);
      tape_stress_yy.insert(tape_stress_yy.begin(), 0);
      tape_stress_zz.insert(tape_stress_zz.begin(), 0);
      tape_stress_yz.insert(tape_stress_yz.begin(), 0);
      tape_stress_xz.insert(tape_stress_xz.begin(), 0);
      tape_stress_xy.insert(tape_stress_xy.begin(), 0);

      // zero!
      vecDispLoad.insert(vecDispLoad.begin(), 0);
      vecTrForce.insert(vecTrForce.begin(), 0);


}
// *****************************************************************************
// *****************************************************************************
void extractdata(mesh mymesh, int body, int avg_box, double vol,
            expression strain_tot, expression stress_corr,
            vector<double> thick_lays, double w, double len)
{
    // this function extracts:
    // 1)strain 2)stress 3)plastic strain 4)non-Mises 5)yield stress 6)whole tape
    // at center point of all layers
    // and store it to the corresponding vector
    // center point of each layer : [w/2, y, len/2]

    double y;         // height of the center point
    double val;       // interpolated value

    int num_lays = thick_lays.size(); // number of layers

    for (int i=0 ; i<num_lays ; i++) {

        y = 0; // resetting y
        // calculating the y
        for (int k=0 ; k<=i ; k++) {
            if (k==i) {
                y = y + thick_lays[k]/2;
            } else {
                y = y + thick_lays[k];
            }
        }

        // strain interpolation @ center and pushing it back
        val = entry(0,0, strain_tot).interpolate(body,   {w/2, y, len/2})[0];
        strain_xx[i]->push_back(val);
        val = entry(1,0, strain_tot).interpolate(body,   {w/2, y, len/2})[0];
        strain_yy[i]->push_back(val);
        val = entry(2,0, strain_tot).interpolate(body,   {w/2, y, len/2})[0];
        strain_zz[i]->push_back(val);
        val = entry(3,0, strain_tot).interpolate(body,   {w/2, y, len/2})[0];
        strain_yz[i]->push_back(val);
        val = entry(4,0, strain_tot).interpolate(body,   {w/2, y, len/2})[0];
        strain_xz[i]->push_back(val);
        val = entry(5,0, strain_tot).interpolate(body,   {w/2, y, len/2})[0];
        strain_xy[i]->push_back(val);

        // stress interpolation @ center and pushing it back
        val = entry(0,0, stress_corr).interpolate(body,    {w/2, y, len/2})[0];
        stress_xx[i]->push_back(val);
        val = entry(1,0, stress_corr).interpolate(body,    {w/2, y, len/2})[0];
        stress_yy[i]->push_back(val);
        val = entry(2,0, stress_corr).interpolate(body,    {w/2, y, len/2})[0];
        stress_zz[i]->push_back(val);
        val = entry(3,0, stress_corr).interpolate(body,    {w/2, y, len/2})[0];
        stress_yz[i]->push_back(val);
        val = entry(4,0, stress_corr).interpolate(body,    {w/2, y, len/2})[0];
        stress_xz[i]->push_back(val);
        val = entry(5,0, stress_corr).interpolate(body,    {w/2, y, len/2})[0];
        stress_xy[i]->push_back(val);

    }

    val = entry(0,0, strain_tot).integrate(avg_box, 4);
    tape_strain_xx.push_back(val/vol);
    val = entry(1,0, strain_tot).integrate(avg_box, 4);
    tape_strain_yy.push_back(val/vol);
    val = entry(2,0, strain_tot).integrate(avg_box, 4);
    tape_strain_zz.push_back(val/vol);
    val = entry(3,0, strain_tot).integrate(avg_box, 4);
    tape_strain_yz.push_back(val/vol);
    val = entry(4,0, strain_tot).integrate(avg_box, 4);
    tape_strain_xz.push_back(val/vol);
    val = entry(5,0, strain_tot).integrate(avg_box, 4);
    tape_strain_xy.push_back(val/vol);

    val = entry(0,0, stress_corr).integrate(avg_box, 4);
    tape_stress_xx.push_back(val/vol);
    val = entry(1,0, stress_corr).integrate(avg_box, 4);
    tape_stress_yy.push_back(val/vol);
    val = entry(2,0, stress_corr).integrate(avg_box, 4);
    tape_stress_zz.push_back(val/vol);
    val = entry(3,0, stress_corr).integrate(avg_box, 4);
    tape_stress_yz.push_back(val/vol);
    val = entry(4,0, stress_corr).integrate(avg_box, 4);
    tape_stress_xz.push_back(val/vol);
    val = entry(5,0, stress_corr).integrate(avg_box, 4);
    tape_stress_xy.push_back(val/vol);

    //
}
// *****************************************************************************
// *****************************************************************************
void savedata(vector<string> tags)
{
    // this function writes all stored data in txt files

    // temperature data
    writevector("./matlab_folder/results/T.txt", vecT, '\n');
    // displacement data
    writevector("./matlab_folder/results/disp.txt", vecDispLoad, '\n');
    // transverse force data
    writevector("./matlab_folder/results/trans_force.txt", vecTrForce, '\n');

    // writing all data in all layers
    for (int i=0; i<tags.size() ; i++) {

        // writing strain data in txt files
        writevector("./matlab_folder/results/" + tags[i] + "_strain_xx.txt", *strain_xx[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_strain_yy.txt", *strain_yy[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_strain_zz.txt", *strain_zz[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_strain_yz.txt", *strain_yz[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_strain_xz.txt", *strain_xz[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_strain_xy.txt", *strain_xy[i], '\n');

        // writing stress data in txt files
        writevector("./matlab_folder/results/" + tags[i] + "_stress_xx.txt", *stress_xx[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_stress_yy.txt", *stress_yy[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_stress_zz.txt", *stress_zz[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_stress_yz.txt", *stress_yz[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_stress_xz.txt", *stress_xz[i], '\n');
        writevector("./matlab_folder/results/" + tags[i] + "_stress_xy.txt", *stress_xy[i], '\n');

    }


    // writing whole tape data in txt files
    writevector("./matlab_folder/results/tape_strain_xx.txt", tape_strain_xx, '\n');
    writevector("./matlab_folder/results/tape_strain_yy.txt", tape_strain_yy, '\n');
    writevector("./matlab_folder/results/tape_strain_zz.txt", tape_strain_zz, '\n');
    writevector("./matlab_folder/results/tape_strain_yz.txt", tape_strain_yz, '\n');
    writevector("./matlab_folder/results/tape_strain_xz.txt", tape_strain_xz, '\n');
    writevector("./matlab_folder/results/tape_strain_xy.txt", tape_strain_xy, '\n');
    //
    writevector("./matlab_folder/results/tape_stress_xx.txt", tape_stress_xx, '\n');
    writevector("./matlab_folder/results/tape_stress_yy.txt", tape_stress_yy, '\n');
    writevector("./matlab_folder/results/tape_stress_zz.txt", tape_stress_zz, '\n');
    writevector("./matlab_folder/results/tape_stress_yz.txt", tape_stress_yz, '\n');
    writevector("./matlab_folder/results/tape_stress_xz.txt", tape_stress_xz, '\n');
    writevector("./matlab_folder/results/tape_stress_xy.txt", tape_stress_xy, '\n');
}
// *****************************************************************************
// *****************************************************************************
