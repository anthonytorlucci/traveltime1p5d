# traveltime1p5d
Travel times and reflection angles for seismic rays in 1.5D models.

### problem definition
Calculate the travel time along raypath given arrays of layer thicknesses $ \left [ q_0, q_1, q_2, ... \right ] $ and velocities $ \left [ v_0, v_1, v_2, ... \right ] $, and distance from the midpoint to source or receiver (referred to here as the half-offset) at the surface $ h $, from the base of the last layer to the surface.

### single layer solution
For a single layer model, the solution for the launch angle is straight forward and calculated deterministically.

$$ \theta_{0} = \tan^{-1}\left ( \frac{h}{q_0} \right ) $$

where $ \theta_{0} $ is the launch angle measured from the vertical for the first layer (`index 0`).

### two-layer solution
However, finding the launch angle from the base of a two layer system cannot be solved deterministically and requires a heuristic since the solution saught $ \theta_1 $ is implicitly defined.

Start with Snell's law for refraction

$$ \frac{\sin\theta_{1} }{v_1} = \frac{\sin\theta_{0} }{v_0} $$

which implies

$$ \theta_{0} = \sin^{-1}\left ( \frac{v_0 \sin\theta_1}{v_1} \right ) $$

where $ \sin\theta_{1} $ is the launch angle from the base of the second layer (`index 1`) and $ \sin\theta_{0} $ is the angle of the refracted ray measured from vertical.

The next step is to define the horizontal distance traveled in each layer.

$$ \zeta_1 = q_1 \tan \left ( \theta_1 \right ) $$
and
$$ \zeta_0 = q_0 \tan \left ( \theta_0 \right ) $$

The final step is to find the angle $ \theta_1 $ that minimizes the error between the user defined half-offset $ h $ and the total horizontal distance traveled $ \zeta_1 + \zeta_0 $. 

$$ minimize(h-(\zeta_1+\zeta_0)) $$

To do so, "guess" an initial launch angle and calculate the horizontal distance traveled. Use a hill climbing algorithm to minimize the objective function above.

### general n-layer solution
Given the incident angle from the base of layer $ n $, determine the angle of the refracted ray in layer $ n-1 $.

$$ \theta_{n-1} = \sin^{-1}\left ( \frac{v_{n-1} \sin\theta_{n}}{v_{n}} \right ) $$

Using horizontal distance traveled in layer $ n $

$$ \zeta_{n} = q_{n} \tan \left ( \theta_{n} \right ) $$

determine the launch angle that minimizes the objective function

$$ minimize(h-(\sum \zeta_{n})) $$

### one-way traveltime
The travel time from the base to the surface along the ray path is

$$ \sum t_{i} = \sum \left ( \frac{r_{i}}{v_{i}} \right ) $$

where $ r_i $ is the length of the ray segment in layer $ i $.

$$ r_{i}=\frac{q_{i}}{\cos\theta_{i}} $$



