# hopping-hoop-simulator
Simulate and animate a hopping hoop! Energy conservation is confirmed throughout.

My only assumptions are (1) flat ground, (2), if possible, rolling-without-sliding initial velocity, and (3) the edge mass is exactly at the radius.

For fun, the code will, by default, continue after landing a hop. I assume that the time of the collision (the time until vy_center is 0) is negligibly small. Sometimes, another hop occurs just after landing!

Hopping hoop...  
[https://youtu.be/ETRpkp03stQ](https://youtu.be/ETRpkp03stQ)

Some papers care about skidding vs slipping, which I guess is found by comparing the direction of the center of the hoop's velocity, *v\_x*, to the direction of the force of kinetic friction. My code prints both directions, but I prefer to just call all of it sliding. Then, if the hoop starts with the edge mass being straight down, there are two simple cases...
* rolling without sliding → slide → rolling without sliding → slide → ...
* rolling without sliding → slide → hop

Run the Python file, hop.py! I use Matplotlib for graphs and animation, and I use SciPy for the numerical solvers and integrators. Change the parameters in the parameter section near the top of the file! Since SciPy's solve_ivp() solves the differential equations, you don't even have to know about differential equations to follow all of the code!

# general math
The best derivation of the equation of motion (the differential equation that solves to give *θ*(*t*)) is to apply Newton's 2nd law to the edge mass, *M*, and the hoop of mass *m* separately, then use Newton's 3rd law to relate the contact force between them. Finally, apply Newton's 2nd law for rotation (net torque = *I* *α*) to the *hoop*. To get the acceleration of *M*, its velocity is -*r* *ω* sin*θ* + *v\_x* in the *x* direction and *r* *ω* cos*θ* in the *y* direction, where *v\_x* = - *r* *ω* is true when rolling without sliding.

Another equivalent derivation method is to apply Newton's 2nd law and Newton's 2nd law for rotation to the *center of mass* (CM). You can always use any CM as its pivot even if it is accelerating. Whenever using Newton's 2nd law, the acceleration is the acceleration of the object's CM. 

Whatever you do, I recommend that you check your rolling-without-sliding equation of motion with what you get from Lagrangian mechanics. To do this, define *L*, the Lagrangian, as kinetic energy minus potential energy, and write is so that *θ* and its time derivative (*ω*) are the only time-dependent variables that appear. Then, ∂*L*/∂*θ* = d(∂*L*/∂*ω*)/dt gives the equation of motion.

Keep in mind that, whenever calculating kinetic energy, it only equals rotational + translational if each are describing the CM. So, for convenience, for everything but the hop, I find the kinetic energy of each object separately then add them together.

During the hop, it is best to use the CM of the combined object for Newton's 2nd law and for calculating energy. You will need the CM formula to get *r*\_CM = *M*/(*m*+*M*) *r*, and you will need the parallel axis theorem to get *I*\_CM = *I* + *M* *r*^2 - (*m*+*M*) *r*\_CM^2, where *I* is the moment of inertia of just the hoop, but *I*\_CM is of the combined object.

# sliding and hopping math
Any sliding starts when *μ*\_*s* *F*_normal - | *F*_friction | becomes negative. The equation of motion is a mess when sliding because terms that would have cancelled no longer cancel. From your system of equations, instead of solving for *F*_friction, you now need to solve for d(*v\_x*)/d*t*.

An interesting part of the calculation is getting the sign of *F*_friction. A way that I could think to do it is try each sign, then use the one and only one that has the sign of *F*_friction be the opposite of the sign of d(*v\_x*)/d*t* + *α* *r*. If *v\_x* + *ω* *r* ever equals 0, either (1) the sign of *F*_friction changes, or (2), if *μ*\_*s* *F*_normal - | *F*_friction | is positive (*F*_friction calculated as if rolling without sliding), the hoop goes back to rolling without sliding. If you were just rolling without sliding, an easier way to find the sign of *F*_friction is that it will not change as sliding starts, so my code could be simplified.

Any hop starts whenever *F*_normal becomes negative, then the CM moves at its previous when-rolling velocity then undergoes projectile motion until the center (not the CM) of the hoop reaches its original height again. The rotational speed, *ω*, does not change when the hop begins. You can use conservation of energy or conservation of angular momentum to prove this, and I highly recommend that, to do this proof, choose the frame of reference that is moving at *v\_x* relative to the ground with origin at the center of mass.

As for landing a hop, there are comments in the code with all of the mathematical details. The basic idea is, since I assume that negligible time elapses, to treat the vertical impulse like time.

The uniquely tricky part of simulating a hopping hoop was realizing that sliding always occurs before a hop, even when the coefficient of static friction, *μ\_s*, is infinite. When first writing the code I initially assumed an infinite *μ\_s* so that I would not have to worry about sliding yet. I reasoned that *μ*\_*s* *F*_normal - | *F*_friction | would never become negative if *μ\_s* is infinite. However, I would get a negative acceleration of the center of the hoop at the very beginning of any hop, which is a contradiction with the hoop hopping. The correct way to take a limit is to have a physical (in this case, non infinite) value, solve the problem, *then* take the limit. When doing this, we can see that, assuming nonzero *F*_friction, *μ*\_*s* *F*_normal - | *F*_friction | will always become negative before *F*_normal becomes negative. Taking the limit as *μ\_s* becomes infinite, sliding and hopping both occur at the same time, but we now see that *sliding* is what results, not hopping.
