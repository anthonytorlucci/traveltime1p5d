"""calculate the travel time for 1.5 D model"""

# standard
from types import FunctionType

# third party
import numpy

def snells_law(theta_i:float, v_i:float, v_r:float):
    """Return the angle for the refracted ray in radians

    Parameters
    ----------
    theta_i: float
        incident angle in radians
    v_i: float
        velocity in incident ray medium
    v_r: float
        velocity in refracted ray medium
    
    Returns
    -------
    float
        angle relative to normal of refracted ray in radians
    """
    return numpy.arcsin(v_r * numpy.sin(theta_i) / v_i)

def horizontal_distance_traveled(theta:float, v:numpy.ndarray, q:numpy.ndarray):
    """Calculate the total horizontal distance traveled for a given launch angle in 1.5D media (assumes flat layer-cake model). 
    Velocity units and layer thickness units for length should be consistent, i.e. meters.

    This is a simple example of the shooting method: a ray is "shot" from a point at the base of the layers from some angle, theta,
    relative to the normal and traversed to the surface. Using Snells law of refraction, the angle of the transmitted ray 
    is calculated for the layer (n-1) from which the horizontal distance is calculated. The total horizontal distance traveled is 
    the sum of the horizontal distance traversed in each layer.

    Parameters
    ----------
    theta: float
        launch angle from the base of the last layer in radians
    v: numpy.ndarray
        array of velocities from surface down (bottom posted)
    q: numpy.ndarray
        array of layer thicknesses from surface down    
    
    Returns
    -------
    float
        distance along surface from midpoint to source or receiver
    """
    # function used in local minmization heuristic
    N = len(v)
    assert(len(q) == N)
    if N == 1:
        h = numpy.array([q[0] * numpy.tan(theta)])
    else:
        phi = numpy.zeros(shape=(N))  # angles at each interface 
        phi[-1] = theta
        h = numpy.zeros(shape=(N))
        #h[-1] = numpy.array([q[-1] * numpy.tan(theta)])
        h[-1] = q[-1] * numpy.tan(theta)
        for n in range(N-1,0,-1):
            phi[n-1] = snells_law(theta_i=phi[n], v_i=v[n], v_r=v[n-1])
            #h[n-1] = numpy.array([q[n-1] * numpy.tan(phi[n-1])])
            h[n-1] = q[n-1] * numpy.tan(phi[n-1])
        
    return numpy.sum(h)

def oneway_traveltime(theta:float, v:numpy.ndarray, q:numpy.ndarray):
    """Calculate the traveltime from the base of the model (v,q) to the surface given the launch angle theta.
    Velocity units and layer thickness units for length should be consistent, i.e. meters. The returned value has the same units as 
    the time units in velocity.

    This is a simple example of the shooting method: a ray is "shot" from a point at the base of the layers from some angle, theta,
    relative to the normal and traversed to the surface. Using Snells law of refraction, the angle of the transmitted ray 
    is calculated for the layer (n-1) from which the travel time is calculated from the length of the ray and the velocity in layer. 
    The total travel time is the sum of the travel time in each layer.

    Parameters
    ----------
    theta: float
        launch angle from the base of the last layer in radians
    v: numpy.ndarray
        array of velocities from surface down
    q: numpy.ndarray
        array of layer thicknesses from surface down
    
    Returns
    -------
    float
        distance along surface from midpoint to source or receiver
    """
    N = len(v)
    assert(len(q) == N)
    if N == 1:
        r = q[0] / numpy.cos(theta)  # length of ray in layer
        t = numpy.array([r / v[0]])  # travel time
    else:
        phi = numpy.zeros(shape=(N))  # angles at each interface 
        phi[-1] = theta
        t = numpy.zeros(shape=(N))
        r = q[-1] / numpy.cos(theta)
        t[-1] = r / v[-1]
        for n in range(N-1,0,-1):
            phi[n-1] = snells_law(theta_i=phi[n], v_i=v[n], v_r=v[n-1])
            r = q[n-1] / numpy.cos(phi[n-1])
            t[n-1] = r / v[n-1]
        
    return numpy.sum(t)

# -----------------------------------------------------------------------------
# https://machinelearningmastery.com/iterated-local-search-from-scratch-in-python/

# converting the two-dimensional domain given in example article above to one-dimensional domain for this travel time problem

# check if a point is within the bounds of the search
def in_bounds1D(point, bounds):
    if point < bounds[0] or point > bounds[1]:
        return False
    return True
 
# hill climbing local search algorithm
def hillclimbing1D(objective, bounds, n_iterations, step_size, start_pt):
    # store the initial point
    solution = start_pt
    # evaluate the initial point
    solution_eval = objective(solution)
    # run the hill climb
    for _ in range(n_iterations):
        # take a step
        candidate = None
        while candidate is None or not in_bounds1D(candidate, bounds):
            candidate = solution + numpy.random.randn() * step_size
        # evaluate candidate point
        candidate_eval = objective(candidate)
        # check if we should keep the new point
        if candidate_eval <= solution_eval:
            # store the new point
            solution, solution_eval = candidate, candidate_eval
    return [solution, solution_eval]

# iterated local search algorithm
def iterated_local_search1D(objective:FunctionType, bounds:list, n_iter:int, step_size:float, n_restarts:int, p_size:float):
    """General Iterated local search in 1D media.

    Parameters
    ----------
    objective: FunctionType
        The objective function to minimize
    bounds: list
        domain of the search space defined as the lower bound and upper bound
    n_iter: int
        total number of iterations
    step_size: float
        maximum step size
    n_restarts: int
        total number of random restarts
    p_size: float
        perturbation step size
    
    Returns
    -------
    list
        best value in the domain and the solution evaluated at best
    """
    # define starting point
    best = None
    while best is None or not in_bounds1D(best, bounds):
        best = bounds[0] + numpy.random.rand() * (bounds[1] - bounds[0])
    # evaluate current best point
    best_eval = objective(best)
    # enumerate restarts
    for n in range(n_restarts):
        # generate an initial point as a perturbed version of the last best
        start_pt = None
        while start_pt is None or not in_bounds1D(start_pt, bounds):
            start_pt = best + numpy.random.randn() * p_size
        # perform a stochastic hill climbing search
        solution, solution_eval = hillclimbing1D(objective, bounds, n_iter, step_size, start_pt)
        # check for new best
        if solution_eval < best_eval:
            best, best_eval = solution, solution_eval
            print('Restart %d, best: f(%s) = %.5f' % (n, best, best_eval))
    return [best, best_eval]

# -----------------------------------------------------------------------------

class TravelTime1D(object):
    def __init__(self, 
        pv:numpy.ndarray,
        dz:numpy.ndarray, 
        h=1000.0, 
        bounds=[0, 0.5*numpy.pi], 
        n_iter=1000,
        step_size=0.01,
        n_restarts=30,
        p_size=1.0
        ):
        #self._df = df
        # TODO assert pv and dz are 1D vectors and have same length or number of elements
        self._pv = pv
        self._dz = dz
        self._h = h  # half-offset for which to calculate the launch angle and travel time
        # variables for luanch angle iterated search algorithm
        self._bounds = bounds  # minimum and maximum values allowed for launch angle; bounds or constraints on the solution space
        self._n_iter = n_iter  # number of iterations ... TODO more information
        self._step_size = step_size  # step size ...TODO more information
        self._n_restarts = n_restarts  # number of restarts for stochastic hill climbing algorithm
        self._p_size = p_size  # perturbation size ...TODO more information
        # results of launch angle search
        best, _ = self._find_launch_angle1d()
        self._best_theta = best  # the best launch angle 
        self._best_h = horizontal_distance_traveled(theta=self._best_theta, v=self._pv, q=self._dz)  # the half offset given for best_theta  

    def _launchangle1d_objective_funtion(self, a:float):
        """The objective funtion is simply the absolute value(so it's always 
        positive) of the difference between the guessed half-offset and the 
        actual half-offset given the launch angle of the ray relative the 
        z-axis perpendicular to the surface.

        Parameters
        ----------
        a: float
            launch angle
        
        Returns
        -------
        float
            squared difference
        """
        #--pv = self._df['pwave_velocity'].values
        #--dz = self._df['layer_thickness'].values
        g = horizontal_distance_traveled(theta=a, v=self._pv, q=self._dz)
        return numpy.abs(self._h - g)

    
    def _find_launch_angle1d(self):
        """Determine the angle of a ray shot from the base of the last layer 
        that minimizes the squared difference between the calculated 
        horizontal distance traveled and the half-offset of the object.
        """ 
        print("starting iterated local search for launch angle ...")
        best, score = iterated_local_search1D(
            objective=self._launchangle1d_objective_funtion, 
            bounds=self._bounds, 
            n_iter=self._n_iter,
            step_size=self._step_size,
            n_restarts=self._n_restarts,
            p_size=self._p_size)
        hbest = horizontal_distance_traveled(theta=best, v=self._pv, q=self._dz)
        print("best solution found with half-offset: {}".format(hbest))
        return best, score

    def traveltime1D(self):
        """Get the travel time from the base of the layer to the surface.

        Returns
        -------
        float
            The one-way travel time from the base of the model to the surface.
        """
        #--pv = self._df['pwave_velocity'].values
        #--dz = self._df['layer_thickness'].values
        # if self._best_theta not None
        return oneway_traveltime(theta=self._best_theta, v=self._pv, q=self._dz)

    def launch_angle(self):
        "Returns the launch angle."
        return self._best_theta