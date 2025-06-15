from __future__ import print_function
import numpy as np
from scipy.optimize import minimize, differential_evolution
import random

pi = np.pi

def metric_tensor(l):
	a,b,c,alpha,beta,gamma = l
	Z = np.array([[a*a, a*b*np.cos(gamma * (pi/180)), a*c*np.cos(beta * (pi/180))],
			      [a*b*np.cos(gamma * (pi/180)), b*b, b*c*np.cos(alpha * (pi/180))],
			  	  [a*c*np.cos(beta * (pi/180)), b*c*np.cos(alpha * (pi/180)), c*c]])
	return Z

def objective(V, ncra, ncca, Alpha, ne, nv, Bstar_inv, SBU_IP):
	
	P = V[0:6]
	Z = metric_tensor(P)
	X = np.array(V[6:])
	X = X.reshape(ncra,ncca)
	var_alpha = np.r_[Alpha[0:ne-nv+1,:], X]

	omp = np.dot(Bstar_inv, var_alpha)
	g = np.dot(omp, np.dot(Z,omp.T))

	O1 = 0
	for sbu in SBU_IP:
		O2 = 0
		for i in sbu[1]:
			O2 = O2 + (i[2] - g[i[0] - 1][i[1] - 1])**2
		O1 = O1 + O2

	return O1


def scale(
	all_SBU_coords,
	a, b, c,
	ang_alpha, ang_beta, ang_gamma, max_le,
	num_vertices, Bstar, alpha, num_edges,
	FIX_UC, SCALING_ITERATIONS, PRE_SCALE, MIN_CELL_LENGTH, OPT_METHOD):
	"""
	Make cell.

	Al parecer es la funcion que inserta en la celda los edges y nodos.

	Parameters:
	-----------
	all_SBU_coords : list of tuples.
		Lista de bloques de construccion donde estan contenidos los nodos. Esto contiene las coordenadas
		de lols atomos dummy X.
	"""
	print("REVISAR: VAMOS A VER ESTA VAINA")
	# compute the max distance between the node center and edges.
	# this goint to be scaled
	max_length = 0
	for line in all_SBU_coords:
		for length in [np.linalg.norm(s[1]) for s in line[1]]:
			if length > max_length:
				max_length = length

	if PRE_SCALE is None:
		scale_guess = 1.0
	else:
		scale_guess = (max_length/max_le) * PRE_SCALE

	all_SBU_ip = []
	all_SBU_ip_append = all_SBU_ip.append
	for sbu in all_SBU_coords:
		SBU_ip = []
		SBU_ip_append = SBU_ip.append
		w = len(sbu[1])
		for i in range(w):
			ivec = sbu[1][i][1]
			iind = sbu[1][i][0]
			for j in range(i, w):
				jvec = sbu[1][j][1]
				jind = sbu[1][j][0]
				dot = np.dot(ivec, jvec)
				SBU_ip_append([iind, jind, dot])
		all_SBU_ip_append((sbu[0], SBU_ip))
	ncra = num_vertices - 1
	ncca = 3

	covars_values = []
	covars_values_append = covars_values.append
	for i in range(ncra):
		for j in range(ncca):
			covars_values_append(0)
			
	# topology cell
	ucvars = [a, b, c, ang_alpha, ang_beta, ang_gamma]

	if np.any(FIX_UC):
		uc_bounds = []
		uc_bounds_append = uc_bounds.append
		for f, p in zip(FIX_UC, ucvars):
			if f:
				uc_bounds_append((p, p))
			else:
				uc_bounds_append((0, None))
		uc_bounds = tuple(uc_bounds)
	else:
		max_a = 2 * scale_guess * a
		max_b = 2 * scale_guess * a
		max_c = 2 * scale_guess * a
		uc_bounds = ((MIN_CELL_LENGTH, max_a), (MIN_CELL_LENGTH, max_b), (MIN_CELL_LENGTH, max_c), (20, 160), (20, 160), (20, 160))

	init_variables = [scale_guess * a, scale_guess * b, scale_guess * c, ang_alpha, ang_beta, ang_gamma] + covars_values
	x_bounds = tuple([(-10.0, 10.0) for x in covars_values])
	bounds = uc_bounds + x_bounds

	Bstar_inv = np.linalg.inv(Bstar)
	# print('scaling unit cell and vertex positions...')

	if OPT_METHOD == 'L-BFGS-B':
		print('optimizing with local minimization algorithm L-BFGS-B...')
	elif OPT_METHOD == 'differential_evolution':
		print('optimizing with global minimization algorithm differential_evolution...')
	else:
		raise ValueError('optimization method', OPT_METHOD, 'is not implemented')

	niter = SCALING_ITERATIONS
	# print("niter:", niter)
	uc_press = 0.05
	covars_perturb = 0.0001
	callbackresults = [[init_variables, objective(
		init_variables, ncra, ncca, alpha, num_edges, num_vertices, Bstar_inv, all_SBU_ip)
	]]

	def callbackF(X):
		var_string = ''
		for val in X:
			var_string += str(val) + '   '
		callbackresults.append([X, objective(X,ncra,ncca,alpha,num_edges,num_vertices,Bstar_inv,all_SBU_ip)])

	for it in range(niter):
		if OPT_METHOD == 'L-BFGS-B':
			res = minimize(
				objective, init_variables,
				args=(ncra, ncca, alpha, num_edges, num_vertices, Bstar_inv, all_SBU_ip),
				method='L-BFGS-B',
				bounds=bounds,
				options={'disp': False},
				callback=callbackF
			)

		elif OPT_METHOD == 'differential_evolution':
			res = differential_evolution(objective, bounds, args=(ncra,ncca,alpha,num_edges,num_vertices,Bstar_inv,all_SBU_ip),
                                    	 polish=True, disp=False, strategy='randtobest1bin')

		uc_params = res.x[0:6]
		uc_params = [i - (i * uc_press) for i in uc_params]
		mult = [random.choice([-1, 1]) for i in range(len(res.x[6:]))]
		covars = [i + j * (i * covars_perturb) for i, j in zip(res.x[6:], mult)]
		init_variables = uc_params + covars

	if niter != 0:	
		sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = res.x[0:6]
		sc_covar = res.x[6:].reshape(ncra,ncca)
	else:
		init_variables = np.asarray(init_variables)
		sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = init_variables[0:6]
		sc_covar = init_variables[6:].reshape(ncra,ncca)

	ab = [a/b, sc_a/sc_b]
	ac = [a/c, sc_a/sc_c]
	bc = [b/c, sc_b/sc_c]
	alpha = [ang_alpha, sc_alpha]
	beta = [ang_beta, sc_beta]
	gamma = [ang_gamma, sc_gamma]
	covar = [np.average(abs(np.array(callbackresults[0][0][6:]) - np.array(callbackresults[-1][0][6:])))]
	final_obj = [callbackresults[-1][1]]

	print('The final objective function value is', np.round(final_obj[0],3))
	print()

	scaling_data = [ab, ac, bc, alpha, beta, gamma, covar, final_obj]

	return(sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma,sc_covar,Bstar_inv,max_length,callbackresults,ncra,ncca,scaling_data)


class UnitCellScaler:
    """
    Class for optimizing unit cell dimensions and vertex positions based on
    building block coordinates.
    """

    def __init__(self,
                 all_SBU_coords, cell_params, max_le,
                 num_vertices, Bstar, alpha, num_edges,
                 FIX_UC, SCALING_ITERATIONS, PRE_SCALE,
                 MIN_CELL_LENGTH, OPT_METHOD):

        self.all_SBU_coords = all_SBU_coords
        self.a = cell_params["aL"]
        self.b = cell_params["bL"]
        self.c = cell_params["cL"]
        self.ang_alpha = cell_params["alpha"]
        self.ang_beta = cell_params["beta"]
        self.ang_gamma = cell_params["gamma"]
        self.max_le = max_le
        self.num_vertices = num_vertices
        self.Bstar = Bstar
        self.alpha = alpha
        self.num_edges = num_edges
        self.FIX_UC = FIX_UC
        self.SCALING_ITERATIONS = SCALING_ITERATIONS
        self.PRE_SCALE = PRE_SCALE
        self.MIN_CELL_LENGTH = MIN_CELL_LENGTH
        self.OPT_METHOD = OPT_METHOD

    def _compute_max_length(self):
        max_length = 0
        for line in self.all_SBU_coords:
            for length in [np.linalg.norm(s[1]) for s in line[1]]:
                if length > max_length:
                    max_length = length
        return max_length

    def _build_inner_products(self):
        all_SBU_ip = []
        for sbu in self.all_SBU_coords:
            SBU_ip = []
            w = len(sbu[1])
            for i in range(w):
                ivec = sbu[1][i][1]
                iind = sbu[1][i][0]
                for j in range(i, w):
                    jvec = sbu[1][j][1]
                    jind = sbu[1][j][0]
                    dot = np.dot(ivec, jvec)
                    SBU_ip.append([iind, jind, dot])
            all_SBU_ip.append((sbu[0], SBU_ip))
        return all_SBU_ip

    def optimize(self):
        print("[INFO] Starting unit cell scaling optimization")

        max_length = self._compute_max_length()
        scale_guess = 1.0 if self.PRE_SCALE is None else (max_length / self.max_le) * self.PRE_SCALE

        all_SBU_ip = self._build_inner_products()

        ncra = self.num_vertices - 1
        ncca = 3
        covars_values = [0] * (ncra * ncca)

        ucvars = [self.a, self.b, self.c, self.ang_alpha, self.ang_beta, self.ang_gamma]
        if np.any(self.FIX_UC):
            uc_bounds = [(p, p) if f else (0, None) for f, p in zip(self.FIX_UC, ucvars)]
        else:
            max_a = 2 * scale_guess * self.a
            max_b = 2 * scale_guess * self.b
            max_c = 2 * scale_guess * self.c
            uc_bounds = [(self.MIN_CELL_LENGTH, max_a),
                         (self.MIN_CELL_LENGTH, max_b),
                         (self.MIN_CELL_LENGTH, max_c),
                         (20, 160), (20, 160), (20, 160)]

        init_variables = [scale_guess * v for v in [self.a, self.b, self.c]] + ucvars[3:] + covars_values
        x_bounds = [(-10.0, 10.0)] * len(covars_values)
        bounds = tuple(uc_bounds + x_bounds)

        Bstar_inv = np.linalg.inv(self.Bstar)

        print(f"[INFO] Optimization method: {self.OPT_METHOD}")
        if self.OPT_METHOD not in ['L-BFGS-B', 'differential_evolution']:
            raise ValueError(f"Optimization method {self.OPT_METHOD} is not implemented")

        callbackresults = [[init_variables, objective(
            init_variables, ncra, ncca, self.alpha, self.num_edges, self.num_vertices, Bstar_inv, all_SBU_ip)]]

        def callbackF(X):
            callbackresults.append([X, objective(X, ncra, ncca, self.alpha, self.num_edges, self.num_vertices, Bstar_inv, all_SBU_ip)])

        for _ in range(self.SCALING_ITERATIONS):
            if self.OPT_METHOD == 'L-BFGS-B':
                res = minimize(
                    objective, init_variables,
                    args=(ncra, ncca, self.alpha, self.num_edges, self.num_vertices, Bstar_inv, all_SBU_ip),
                    method='L-BFGS-B', bounds=bounds, callback=callbackF, options={'disp': False}
                )
            elif self.OPT_METHOD == 'differential_evolution':
                res = differential_evolution(
                    objective, bounds,
                    args=(ncra, ncca, self.alpha, self.num_edges, self.num_vertices, Bstar_inv, all_SBU_ip),
                    polish=True, disp=False
                )

            uc_params = [v - (v * 0.05) for v in res.x[:6]]
            mult = [random.choice([-1, 1]) for _ in res.x[6:]]
            covars = [v + j * (v * 0.0001) for v, j in zip(res.x[6:], mult)]
            init_variables = uc_params + covars

        final_x = res.x if self.SCALING_ITERATIONS != 0 else np.asarray(init_variables)

        sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma = final_x[:6]
        sc_covar = final_x[6:].reshape(ncra, ncca)

        scaling_data = [
            [self.a / self.b, sc_a / sc_b],
            [self.a / self.c, sc_a / sc_c],
            [self.b / self.c, sc_b / sc_c],
            [self.ang_alpha, sc_alpha],
            [self.ang_beta, sc_beta],
            [self.ang_gamma, sc_gamma],
            [np.average(abs(np.array(callbackresults[0][0][6:]) - np.array(callbackresults[-1][0][6:])))],
            [callbackresults[-1][1]]
        ]

        print(f"[RESULT] Final objective value: {np.round(callbackresults[-1][1], 3)}")

        self.scaled_cell_param = {
            "aL": sc_a,
            "bL": sc_b,
            "cL": sc_c,
            "alpha": sc_alpha,
            "beta": sc_beta,
            "gamma": sc_gamma
        }

        return (self.scaled_cell_param, sc_covar, Bstar_inv, max_length, callbackresults, ncra, ncca, scaling_data)