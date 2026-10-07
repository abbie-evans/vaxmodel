import numpy as np
from scipy.integrate import solve_ivp
from numba import jit

class Model():
    """ 
    Model for the antibody dynamics following vaccination.
    """

    def __init__(self):
        """
        Initialises the model.
        """
        super(Model, self).__init__()

    def simulate(self, times, _params, fit_times=None, doses=1, outputs=3, naive=True, second_dose=None):
        """
        Computes the values in each compartment of the ODE system.

        Parameters
        ----------
        times : list
            List of time points at which vaccination occurs.
        _params : dict
            Dictionary of the model parameters.
        fit_times : list, optional
            List of time points at which to fit the model. If None, uses the
            times provided.
        doses : int
            Number of doses.
        outputs : int
            Number of outputs. If 3, returns total Abs and Ag with time; if 5, returns
            also specific IgG to each epitope.
        naive : bool
            If True, simulates a naive individual; if False, simulates a
            recovered individual with pre-existing immunity.
        second_dose : bool
            If True, returns the antibody levels after the second dose only; if False, returns the antibody levels after n doses.

        Returns
        -------
        list
            Solution of the ODE system at the time points provided.

        """
        dAg = _params['dAg']
        dIC = _params['dIC']
        dAb = _params['dAb']
        s = _params['s']
        phi = _params['phi']
        pp = _params['pp']
        kslow = _params['kslow']
        kA = _params['kA']
        dB = _params['dB']
        nn = _params['nn']
        r = _params['r']
        CC = _params['CC']
        c = _params['c']#k*CC
        dS = _params['dS']
        dL = _params['dL']
        M = doses
        B_0 = _params['B_0']
        mu_m = _params['mu_m']
        m = _params['m']
        mu_S = _params['mu_S']
        mu_L = _params['mu_L']
        r_AM = _params['r_AM']
        dAg_depot = _params['dAg_depot']
        d_n = _params['d_n']
        d_m = _params['d_m']
        dose_p = _params['dose_p']
        dose_o = _params['dose_o']
        Ab_0 = _params['Ab_0']
        MW = _params['MW_IgG']
        dose = _params['dose']

        Ab_ref = 6.022*10**23*10**-6/(MW)
        B_ref = 0.5*10**B_0
        pp_ref = pp*B_ref/Ab_ref
        Ag_ref = 6.022*10**23*10**-12/78300

        dose_ = [((dose/(1+0.01*dose))*10**dose_p)/Ag_ref for i in range(M)]

        params = np.array([dAg, dIC, dAb, s, phi, pp, kA, dB, nn, r, CC, c, dS, dL,
                           B_0, mu_m, m, mu_S, mu_L, r_AM, dAg_depot, d_n, d_m, Ab_0, MW], dtype=np.float64)

        # Initial conditions
        pre_immunity = [_params['B_0'], _params['B_0'], _params['Ab_0'], _params['Ab_0']]
        bx, bs, ax, as_ = pre_immunity

        if len(_params['IC_adm']) == 0:
            IC_adm = [0, 0, 0, 0]
        else:
            IC_adm = _params['IC_adm']

        pp = _params['pp']
        d_Ab = _params['d_Ab']
        lx = d_Ab*ax/pp_ref
        ls = d_Ab*as_/pp_ref
        k_init = 9

        if outputs!=5:
            init = [0, *IC_adm, 1, 1, 0, 0, lx, ls, ax, as_, dose_[0], 0, 0, 0, 0, k_init]
        else:
            init = [0, *IC_adm, 1, 1, 0, 0, lx, ls, ax*15, as_*35, dose_[0], 0, 0, 0, 0, k_init] # Accounting for some prior immunity
            # params[4] = params[4]*2

        if not naive:
            init = [6.28992747e-02, 7.32437311e-02, 5.86279753e-01, 5.86279753e-01,
                    1.78446454e+00, 9.71656375e-01, 9.71656375e-01, 9.84870484e-03,
                    9.84870484e-03, 4.06064136e-04, 4.06064136e-04, 9.98378253e+00,
                    9.98378253e+00, 7.04180389e-05+dose_[0], 3.12541153e-03, 3.12541153e-03,
                    4.10356520e-03, 4.10356520e-03, 9.62399778e+00] # Uses the state following the first dose as the initial condition

        y = []
        y2 = []
        t = []
        ax0 = []
        as0 = []
        sec_dose = []

        # Solve the system of ODEs
        if fit_times is not None:
            if len(fit_times) > 1:
                sol = solve_ivp(
                    lambda t, y: self._right_hand_side(
                        t, y, params),
                    [fit_times[0][0], fit_times[0][-1]], init, t_eval=fit_times[0], method='LSODA')
            else:
                sol = solve_ivp(
                    lambda t, y: self._right_hand_side(
                        t, y, params),
                    [fit_times[0], fit_times[-1]], init, t_eval=fit_times, method='LSODA')
        else:
            sol = solve_ivp(
                lambda t, y: self._right_hand_side(
                    t, y, params),
                [times[0], times[1]], init, t_eval=np.linspace(times[0], times[1], 366), method='LSODA')

        y.append(sol.y[11]+sol.y[12])
        y2.append(sol.y[0])
        t.append(sol.t)
        if outputs==5:
            ax0.append(sol.y[11])
            as0.append(sol.y[12])

        # Boosts
        mm=2 # This is to deal with M=1
        if _params['Het']==0:
            while(mm<=M): # Does not execute boosts if M=1
                init = sol.y[:,-1]
                init[1:5] = [x + y for x, y in zip(init[1:5], IC_adm)]
                init[13] = init[13] + dose_[mm-1]
                if fit_times is not None:
                    if len(fit_times) > 1:
                        sol = solve_ivp(
                            lambda t, y: self._right_hand_side(
                                t, y, params),
                            [fit_times[mm-1][0], fit_times[mm-1][-1]], init, t_eval=fit_times[mm-1], method='LSODA')
                    else:
                        sol = solve_ivp(
                            lambda t, y: self._right_hand_side(
                                t, y, params),
                            [fit_times[0], fit_times[-1]], init, t_eval=fit_times, method='LSODA')
                else:
                    sol = solve_ivp(
                        lambda t, y: self._right_hand_side(
                            t, y, params),
                        [times[mm-1], times[mm]], init, t_eval=np.linspace(times[mm-1], times[mm], 366), method='LSODA')
                y.append(sol.y[11]+sol.y[12])
                y2.append(sol.y[0])
                t.append(sol.t)
                sec_dose.append(sol.y[11]+sol.y[12])
                if outputs==5:
                    ax0.append(sol.y[11])
                    as0.append(sol.y[12])
                mm = mm+1

        t = np.array([item for sublist in t for item in sublist])
        ab = np.array([item for sublist in y for item in sublist])
        ant = np.array([item for sublist in y2 for item in sublist])
        hx0 = np.array([item for sublist in ax0 for item in sublist])
        h0s = np.array([item for sublist in as0 for item in sublist])
        sec_dose = np.array([item for sublist in sec_dose for item in sublist])
        if outputs!=5 and second_dose is None:
            return t, ab, ant
        elif outputs==5 and second_dose is None:
            return t, ab, ant, hx0, h0s
        else:
            return sec_dose

    @staticmethod
    @jit(nopython=True)
    def _right_hand_side(t, y, params):
        """
        Constructs the RHS of the equations of the system of ODEs for given a
        time point.

        Parameters
        ----------
        t : float
            Time point at which we compute the evaluation.
        y : numpy.array
            Array of all the compartments of the ODE system.
        params : numpy.array
            Array of all the model parameters.

        Returns
        -------
        numpy.array
            Matrix representation of the RHS of the ODEs system.

        """
        # Split compartments into their types

        # Free antigen, Complement-bound ICs
        H_XS, C_XS = y[0], y[1]
        # Antibody bound ICs
        H_OS, H_XO, H_OO = y[2], y[3], y[4]
        # Resting naive B cells
        B_X, B_S = y[5], y[6]
        # Short-lived plasma cells
        S_X, S_S = y[7], y[8]
        # Long-lived plasma cells
        L_X, L_S = y[9], y[10]
        # Antibodies
        A_X, A_S = y[11], y[12]
        # Depot
        A = y[13]
        # Memory B cells
        M_X, M_S = y[14], y[15]
        # Activated naive B cells
        N_X, N_S = y[16], y[17]
        # Affinity
        k = y[18]

        # Read parameters
        dAg, dIC, dAb, s, phi, pp, kA,\
        dB, nn, r, CC, c, dS, dL, B_0, mu_m, m, mu_S, mu_L, r_AM, dAg_depot, d_n, d_m, Ab_0, MW = params

        dAg = 10**dAg
        kA = 10**kA
        mu_L = 10**mu_L
        mu_S = 10**mu_S
        dS = 10**dS
        
        Ag_B_X  = C_XS + H_XO
        Ag_B_S  = C_XS + H_OS
        k_max = 11
        Ab_ref = (6.022*10**23*10**-6)/MW
        B_ref = 0.5*10**B_0
        pp_ref = pp*B_ref/Ab_ref
        Ag_ref = (6.022*10**23*10**-12)/78300

        dydt = np.array((
            kA*A - dAg*H_XS - 10**(k-np.log10(6.022*10**23))*10**CC*H_XS - (10**(k-np.log10(6.022*10**23)))*Ab_ref*H_XS*(A_X + A_S),
            10**(k-np.log10(6.022*10**23))*10**CC*H_XS - 10**(k-np.log10(6.022*10**23))*Ab_ref*C_XS*(A_X + A_S) - dAg*C_XS,
            10**(k-np.log10(6.022*10**23))*Ab_ref*((H_XS+C_XS)*A_X - H_OS*A_S) - dIC*H_OS,
            10**(k-np.log10(6.022*10**23))*Ab_ref*((H_XS+C_XS)*A_S - H_XO*A_X) - dIC*H_XO,
            10**(k-np.log10(6.022*10**23))*Ab_ref*(H_OS*A_S + H_XO*A_X) - dIC*H_OO,
            dB - (s*B_X*Ag_B_X**nn)/(((10**(phi-k+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_X**nn) - dB*B_X,
            dB - (s*B_S*Ag_B_S**nn)/(((10**(phi-k+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_S**nn) - dB*B_S,
            mu_S*N_X + s*M_X*Ag_B_X**nn/(((10**(phi-k-m+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_X**nn) - (dS+mu_L)*S_X,
            mu_S*N_S + s*M_S*Ag_B_S**nn/(((10**(phi-k-m+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_S**nn) - (dS+mu_L)*S_S,
            mu_L*S_X - dL*L_X,
            mu_L*S_S - dL*L_S,
            pp_ref*(S_X + L_X) - dAb*A_X,
            pp_ref*(S_S + L_S) - dAb*A_S,
            -kA*A - dAg_depot*A,
            mu_m*N_X - d_m*M_X - (s*M_X*Ag_B_X**nn)/(((10**(phi-k-m+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_X**nn),
            mu_m*N_S - d_m*M_S - (s*M_S*Ag_B_S**nn)/(((10**(phi-k-m+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_S**nn),
            s*B_X*Ag_B_X**nn/(((10**(phi-k+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_X**nn) - (d_n + mu_S + mu_m)*N_X,
            s*B_S*Ag_B_S**nn/(((10**(phi-k+np.log10(6.022*10**23)))/Ag_ref)**nn + Ag_B_S**nn) - (d_n + mu_S + mu_m)*N_S,
            r_AM*k*(1-k/k_max), # Affinity maturation
            ))

        return dydt
