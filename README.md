# hopping-hoop-simulator
Simulate a hopping hoop! My only assumptions are: flat ground and, if possible, rolling-without-slipping initial velocity. Energy conservation is confirmed throughout.

Hopping hoop...  
[https://youtu.be/ETRpkp03stQ](https://youtu.be/ETRpkp03stQ)

Run the Python file, hop.py! I use Matplotlib for graphs and animation, and I use SciPy for the numerical solvers and integrators. Change the parameters in the parameter section near the top of the file!

# general math
The best derivation of the equation of motion (the differential equation that solves to give *θ*(*t*)) is to apply Newton's 2nd law to the edge mass, *M*, and the hoop of mass *m* separately, then use Newton's 3rd law to relate the contact force between them. Finally, apply Newton's 2nd law for rotation (net torque = *I* *α*) to the *hoop*. To get the acceleration of *M*, its velocity is -*r* *ω* sin*θ* + *v\_x* in the *x* direction and *r* *ω* cos*θ* in the *y* direction, where *v\_x* = - *r* *ω* is true when rolling without slipping.

Another equivalent derivation method is to apply Newton's 2nd law and Newton's 2nd law for rotation to the *center of mass* (CM). You can always use the CM as the pivot even if it is accelerating. Whenever using Newton's 2nd law, the acceleration is the acceleration of the object's CM. 

Whatever you do, I highly recommend that you check your rolling-without-slipping equation of motion with what you get from Lagrangian mechanics. Define *L*, the Lagrangian, as kinetic energy minus potential energy, and write is so that *θ* and its time derivative (*ω*) are the only time-dependent variables that appear. Then, ∂*L*/∂*θ* = d(∂*L*/∂*ω*)/dt gives the equation of motion.

Keep in mind that, whenever calculating kinetic energy, it only equals rotational + translational if each are describing the CM. So, for convenience, for everything but the hop, I find the kinetic energy of each object separately then add them together.

During the hop, it is best to use the CM of the combined object for Newton's 2nd law and for calculating energy. You will need the CM formula to get *r*\_CM = *M*/(*m*+*M*) *r*, and you will need the parallel axis theorem to get *I*\_CM = *I* + *M* *r*^2 - (*m*+*M*) *r*\_CM^2, where *I* is the moment of inertia of just the hoop, but *I*\_CM is of the combined object.

# sliding and hopping math
Any sliding starts when *μ*\_*s* *F*_normal - | *F*_friction | becomes negative. The equation of motion is a mess when sliding because terms that would have cancelled no longer cancel. From your system of equations, instead of solving for *F*_friction, you now need to solve for d(*v\_x*)/d*t*. An awkward part of the calculation is getting the sign of *F*_friction. The only way that I could think to do it is try each sign, then use the one and only one that has the sign of *F*_friction be the opposite of the sign of d(*v\_x*)/d*t* + *α* *r*. If *v\_x* + *ω* *r* ever equals 0, either (1), if *μ*\_*s* *F*_normal - | *F*_friction | is positive, the hoop goes back to rolling without slipping, or (2) the sign of *F*_friction changes.

Any hop starts whenever *F*_normal becomes negative, then the CM moves at its previous when-rolling velocity then undergoes projectile motion until the center (not the CM) of the hoop reaches its original height again. The rotational speed, *ω*, does not change when the hop begins. You can use conservation of energy or conservation of angular momentum to prove this, and I highly recommend that, to do this proof, choose the frame of reference that is moving at *v\_x* relative to the ground with origin at the center of mass.


