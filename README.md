# hopping-hoop-simulator
Simulate a hopping hoop! Only assume flat ground and, if possible, rolling-without-slipping initial velocity. Energy conservation is confirmed throughout.

Hopping hoop...  
[https://www.youtube.com/watch?v=ETRpkp03stQ](https://www.youtube.com/watch?v=ETRpkp03stQ)

# math
The best derivation of the equation of motion (the differential equation that solves to give *θ*(*t*)) is to apply Newton's 2nd law to the edge mass, *M*, and the disk of mass *m* separately, then use Newton's 3rd law to relate the contact force between them. Finally, apply Newton's 2nd law for rotation to the *disk*. To get the acceleration of *M*, its velocity is -*r* *ω* sin*θ* + *v\_x* in the *x* direction and *r* *ω* cos*θ* in the *y* direction, where *v\_x* = - *r* *ω* is true when rolling without slipping.

Another good method is to apply Newton's 2nd law and Newton's 2nd law for rotation to the *center of mass* (CM). For an accelerating pivot, you can always use the CM. For Newton's 2nd law, the acceleration is the acceleration of the object's CM. 

Any sliding starts when *μ*\_*s* *F*_normal - | *F*_friction | becomes negative. The equation of motion is a mess when sliding because terms that would have cancelled no longer cancel. From your system of equations, instead of solving for *F*_friction, you now need to solve for d(*v\_x*)/d*t*. An awkward part of the calculation is getting the sign of *F*_friction. The only way that I could do it is try each sign, then use the one and only one that has the sign of *F*_friction be the opposite of *α* + 1/*r* d(*v\_x*)/d*t*.

Any hop starts when *F*_normal becomes negative, then the CM moves at its previous rolling velocity then undergoes projectile motion until the center (not the CM) of the hoop reaches its original height again. You will need the CM formula to get *r*\_CM, and you will need the parallel axis theorem. The rotational speed, *ω*, does not change when the hop begins. You can use conservation of energy or conservation of angular momentum to prove this, and I highly recommend that, to do this proof, choose the frame of reference that is moving at *v\_x* relative to the ground with origin at the center of mass.

Keep in mind that, when calculating kinetic energy, it is only rotational + translational if you do rotational\_CM + translational\_CM. So, for convenience, I find the kinetic energy of each object separately (except during the hop).
