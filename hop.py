#!/usr/bin/env python3.11
#
# Simulate the hopping hoop!
#   https://youtu.be/ETRpkp03stQ
# Assume flat ground but not any initial condition or radial distribution of mass of "hoop".
# If possible, assume that rolling without sliding is the initial condition.
# Assume that the edge mass, M, is exactly at r.
#
# All units are MKS.
#
# Radians are used with counterclockwise being the positive direction.
# theta points to the mass on edge of hoop from the center of the hoop.
# theta=0 points to the right.
#
# I believe that sliding is at least sometimes necessary before a hop for the math to be consistent!
# That is, not checking the rolling-without-sliding ODE for sliding doesn't work.
# That is, even with a huge mu_s, the condition to slide is usually (always?) met before the hop.
#
# The following cases have been observed to occur from this code (starting at -90 degrees)...
#  • rolling without sliding -> slide -> rolling without sliding -> slide
#  • rolling without sliding -> slide -> hop
# Note that the slide phase can have the direction of friction switch one or more times.
# I need to play around with the parameter space some more to see what else can happen!
#
# The logic of my land() function is a bit dubious.
# I included land() because I think it usually works, unless the normal force is negative immediately
#   after vy_center=0 is achieved, which the code checks for.
#
# (c) 2023 Bradley Knockel


from numpy import sin, cos, sqrt, radians, degrees, sign, absolute, diff, linspace, concatenate, isclose, arange, interp
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib import animation

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html



#############################
# Set parameters
#############################
# all must be positive except for...
#  • theta0 and omega0 can be anything
#  • rotational inertia of hoop, I, can be 0 or positive
#  • mu_k can be 0 or positive, but, if it is 0 and a hop doesn't occur, the code will never end!
#  • m can be 0 or positive, but, if m=0, issues arise and stay away from theta being straight down
#  • g can be anything, but, if g <= 0, the calculation of when hop is finished obviously breaks


g = 9.81

# hoop
r = 0.5
m = 0.5
I = 0.9 * m * r * r

# hoop and ground friction coefficients
mu_s = 0.2
mu_k = 0.1

# mass on edge of hoop (kg)
M = 1.0

# initial conditions (rad and rad/s)
omega0 = -11
# Straight down (-90 degrees) and straight up (90 degrees) are great because they
#   can always start without sliding, and make the sign of omega0 trivial.
# Straight down can always begin without hopping and avoids unnecessary sliding.
# When m is very small, straight up has a more meaningful omega0 than straight down does.
theta0 = radians(-90)



# absolute tolerance of simulation
abs_tol = 1E-8

# to slow down the animation compared to real time
delayMultiplier = 2   # a positive integer

## If you search this file for "feel free",
## you will find some other interesting parameters
## that you might want to change.



if mu_k >= mu_s:
  print("  ERROR: your mu_k is not less than your mu_s.")
  exit()   # if mu_k == mu_s, my code to find signFrict doesn't work

if I > m*r*r:
  print("  ERROR: your I is larger than m r^2.")
  exit()



## pre calculate (these variables can be used for sliding and hopping too)
g_over_r = g/r
aTerm = I / (M * r * r) + m/M + 2.0   # this is a constant term in the ODE
weight = (m+M)*g
M_times_r = M*r
m_times_r = m*r
# for center of mass (CM)
r_cm = M*r / (m+M)
I_cm = I + M*r*r - (m+M)*r_cm*r_cm



#############################
# graphing and animation code
#############################

_,ax = plt.subplots(2,2)


def makePlots():

  plt.figure()
  plt.plot( big[0], (degrees(big[1]) + 270)%360 - 270, 'k' )  # I want -90 degrees to be the middle of range
  plt.title("θ (in degrees)")
  plt.figure()
  plt.plot( big[0], big[2], 'k' )
  plt.title("ω")
  plt.figure()
  plt.plot( big[0], big[3], 'k' )
  plt.title("α")
  plt.figure()
  plt.plot( big[0], big[4], 'k' )
  plt.title("v_x")
  plt.figure()
  plt.plot( big[0], cumulative_trapezoid( big[4], x = big[0], initial = 0 ), 'k' )
  plt.title("x")
  plt.figure()
  plt.plot( big[0], big[4] + r*big[2], 'k' )
  plt.title("v_bottom = v_x + r ω")
  if len(big) == 6:   # if hop happened
    plt.figure()
    plt.plot( big[0], big[5], 'k' )
    plt.title("y")

  # the following pauses code execution, so call makePlots() after doing useful print()s
  plt.show()



def animate():

    x_pos = cumulative_trapezoid( big[4], x = big[0], initial = 0 )


    ## get values at more consistent times

    omegaMax = sqrt( ( energyFunc_RollingWithoutSliding((theta0, omega0)) + M*g*r ) / (0.5*I + 0.5*m*r*r) )
    tStep_ms = round(100.0 / omegaMax)   # feel free to change numerator
    if not tStep_ms:
      tStep_ms = 1
    tStep = tStep_ms / 1000.0

    tNew = arange(0.0, big[0][-1], tStep)
    thetaNew = interp( tNew, big[0], big[1] )
    xNew = interp( tNew, big[0], x_pos )
    if len(big) == 6:   # if hop happened
      yNew = interp( tNew, big[0], big[5] )
    else:
      yNew = 0.0*tNew


    coordsList = [(xNew[i],yNew[i],xNew[i]+r*cos(thetaNew[i]),yNew[i]+r*sin(thetaNew[i])) for i in range(len(tNew))]


    # setup figure
    pad = 0.5         # prevents edge of hoop from sometimes being chopped off
    maxInches = 8     # feel free to change!
    xRange = [ min(xNew) - (1+pad)*r, max(xNew) + (1+pad)*r ]
    yRange = [ -(1+pad)*r, max(yNew) + (1+pad)*r ]
    xDelta = xRange[1] - xRange[0]
    yDelta = yRange[1] - yRange[0]
    scale = maxInches / max(( xDelta, yDelta ))
    fig = plt.figure( figsize=(scale*xDelta, scale*yDelta) )
    board = plt.axes(xlim=xRange, ylim=yRange, frameon=False)
    board.axis('off')
    circle = plt.Circle((0.0, 0.0), r, fill=False, color='b')
    arrow, = board.plot([0.0, r*cos(thetaNew[0])], [0.0, r*sin(thetaNew[0])], 'b' ) 
    board.add_patch(circle)


    def animateFrame(coords):
        circle.set( center = (coords[0], coords[1]) )
        arrow.set_data([coords[0], coords[2]], [coords[1], coords[3]])
        return [circle, arrow]



    anim = animation.FuncAnimation(fig, animateFrame, frames=coordsList, interval = tStep_ms*delayMultiplier, repeat_delay=1000, blit=True)   # repeat_delay not working on macOS for some reason



    print("  animation time step in ms =", tStep_ms,"*", delayMultiplier)
    print()

    plt.show()

    ## The following first required `sudo port install ffmpeg` on macOS
    ## On Windows: https://www.geeksforgeeks.org/how-to-install-ffmpeg-on-windows/
    ##   though you only need the ffmpeg.exe file
    #anim.save('animation.mp4', writer = animation.FFMpegWriter(fps=30))     # feel free to uncomment this line!



#############################
# rolling without sliding
#############################



## For the following functions,
##   y[0] is theta; y[1] is omega, the time derivative of theta

def derivatives(t, y):
  return ( y[1], - (y[1]*y[1] + g_over_r) * cos(y[0]) / (aTerm + 2*sin(y[0])) )

def normalForce(t, y):
  return M_times_r * ( derivatives(t,y)[1] * cos(y[0]) - y[1]*y[1]*sin(y[0]) ) + weight

def staticFriction(y):
  alpha = derivatives(0,y)[1]
  return - M_times_r * ( y[1]*y[1]*cos(y[0]) + alpha*(1+sin(y[0])) ) - m_times_r * alpha

def maxFriction_minus_requiredFriction(t, y):   # static friction; magnitudes
  return mu_s * normalForce(t,y) - abs( staticFriction(y) )

normalForce.terminal = True
maxFriction_minus_requiredFriction.terminal = True




def energyFunc_RollingWithoutSliding(y):   # y is [ theta, omega, ... ]
  return (M * r * r * (1 + sin(y[0])) + 0.5 * I + 0.5 * m * r * r) * y[1] * y[1] + M * g * r * sin(y[0])



def rollWithoutSliding(finalT, finalY):


    ## solve the ODE!

    eventTuple = (normalForce, maxFriction_minus_requiredFriction)

    tStop = 10/abs(finalY[1])   # feel free to change!
    tList = linspace(finalT, finalT + tStop, num=101)   # nice for making plot; feel free to change num

    sol = solve_ivp(derivatives, (finalT, finalT + tStop), finalY, method = 'Radau', atol = abs_tol, events = eventTuple, t_eval = tList)



    length = len(sol.t)
    theta = sol.y[0] % radians(360)
    alpha = [derivatives(sol.t[i], sol.y[:,i])[1] for i in range(length)]
    Fn = [normalForce(sol.t[i], sol.y[:,i]) for i in range(length)]
    maxMinusFriction = [maxFriction_minus_requiredFriction(sol.t[i], sol.y[:,i]) for i in range(length)]


    # update the graph
    ax[0,0].plot(sol.t, theta, 'k-', sol.t, sol.y[1], 'r--', sol.t, alpha, 'bs', sol.t, Fn, 'g^', sol.t, maxMinusFriction, 'm+')
    ax[0,0].set(xlabel="t")
    ax[0,0].legend(['θ mod 2π','ω','α','F_n','max - |frict|'])



    # start growing the overall lists if you want to simulate things or graph things later
    big[0] = concatenate(( big[0], sol.t ))    # t
    big[1] = concatenate(( big[1], sol.y[0] )) # theta
    big[2] = concatenate(( big[2], sol.y[1] )) # omega
    big[3] = concatenate(( big[3], alpha ))    # alpha
    big[4] = concatenate(( big[4], -r*sol.y[1] ))  # v_x



    finalT_hop = False
    if sol.t_events[0].size:
        finalT_hop = sol.t_events[0][0]
        finalY_hop = sol.y_events[0][0]

    finalT_slide = False
    if sol.t_events[1].size:
        finalT_slide = sol.t_events[1][0]
        finalY_slide = sol.y_events[1][0]

    # sort out if you hop, slide, or neither
    slide = False
    if finalT_slide and finalT_hop:
        if finalT_slide < finalT_hop:   # what if they are equal?
          finalT = finalT_slide
          finalY = finalY_slide
          slide = True
        else:
          finalT = finalT_hop
          finalY = finalY_hop
    elif finalT_slide:
        finalT = finalT_slide
        finalY = finalY_slide
        slide = True
    elif finalT_hop:
        finalT = finalT_hop
        finalY = finalY_hop
    else:   # neither
        print("  Nothing interesting happened before tStop.")
        print("  sum of |Δtheta| =", degrees( sum(absolute(diff(sol.y[0]))) ), "degrees")
        print("  energy needed to make it over the top is", M*g*r)
        print("  energy =", energyFunc_RollingWithoutSliding(sol.y[:,0]) )
        print()
        makePlots()
        animate()
        exit()


    print("  final t =", finalT)
    print("  final θ =", str(finalY[0]) + ", so " + str(degrees(finalY[0])%360) + "°")
    print("  final ω =", finalY[1])
    print("  final α =", derivatives(finalT, finalY)[1])
    print("  final F_n =", normalForce(finalT, finalY))
    print("  final μ_s F_n - |F_fr| =", maxFriction_minus_requiredFriction(finalT, finalY))
    print("  final F_fr =", staticFriction(finalY) )
    print("  energies:", energyFunc_RollingWithoutSliding(sol.y[:,0]), energyFunc_RollingWithoutSliding(finalY))
    if slide:
      print("  Slide!")
    else:
      print("  Hop!")
    print()


    # add on v_x
    finalY = concatenate( (finalY,  [ - r * finalY[1] ] ) )

    return (finalT, finalY, slide)



#############################
# sliding
#############################


def bottomAcceleration(signFrict, y):
  mu = signFrict * abs(mu_k)
  ratio = (m+M)/M
  ratio2 = mu/ratio
  constTerm = mu * ratio * g_over_r
  constTerm2 = 1 + I/(M*r*r)
  s = sin(y[0])
  c = cos(y[0])
  temp = y[1]*y[1]*s
  alpha = (-mu*temp + constTerm - g_over_r*c + mu*g_over_r*s - ratio2*temp*s + temp*c/ratio) / (constTerm2 - mu*c - ratio2*s*c - s*s/ratio)
  dvxdt_over_r = mu*g_over_r + ratio2*(-temp + alpha*c) + (temp*c/s + alpha*s)/ratio
  return r*(alpha + dvxdt_over_r)


# This is the mechanical energy; y is [ theta, omega, v_x, ... ]
def energyFunc_notHopping(y):
  return 0.5 * (I + M * r * r * cos(y[0])**2) * y[1]**2 + 0.5 * m * y[2]**2 + 0.5 * M * (-r*y[1]*sin(y[0]) + y[2])**2 + M * g * r * sin(y[0])


def rollWithSliding(finalT, finalY):



    ## y is [ theta, omega, v_x, W ], where W doesn't affect any derivatives, and v_x only affects the derivative of W
    ## The only reason to integrate to get W is to check energy conservation. dW/dt = power = |v_bottom * F_fr|

    ## The following 3 functions must be nested functions so that they are in the same scope as...
    ##   signFrict, mu, ratio, ratio2, constTerm, constTerm2

    def derivativesSliding(t, y):
      s = sin(y[0])
      c = cos(y[0])
      temp = y[1]*y[1]*s
      alpha = (-mu*temp + constTerm - g_over_r*c + mu*g_over_r*s - ratio2*temp*s + temp*c/ratio) / (constTerm2 - mu*c - ratio2*s*c - s*s/ratio)
      return ( y[1], alpha, r*(mu*g_over_r + ratio2*(-temp + alpha*c) + (temp*c/s + alpha*s)/ratio), -(r*y[1]+y[2]) * mu * ( M_times_r * (alpha * c - temp) + weight ) )

    def normalForceSliding(t, y):
      #return derivativesSliding(t,y)[3] / ((r*y[1]+y[2]) * mu)    # can divide by 0
      s = sin(y[0])
      c = cos(y[0])
      temp = y[1]*y[1]*s
      alpha = (-mu*temp + constTerm - g_over_r*c + mu*g_over_r*s - ratio2*temp*s + temp*c/ratio) / (constTerm2 - mu*c - ratio2*s*c - s*s/ratio)
      return M_times_r * (alpha * c - temp) + weight

    def speedDiff(t, y):
      return signFrict * (r*y[1] + y[2])   # I added signFrict to always make this negative
    speedDiff.direction = 1                # to prevent the event from always being triggered immediately (this works only when speedDiff is negative)

    normalForceSliding.terminal = True
    speedDiff.terminal = True



    # Get direction of kinetic friction based on sign of initial bottomAcceleration once sliding.
    # signFrict and bottomAcceleration should have opposite signs
    tempPos =  sign(bottomAcceleration(+1.0, finalY)) == -1    # try  signFrict = 1
    tempNeg =  sign(bottomAcceleration(-1.0, finalY)) == 1     # try  signFrict = -1
    if tempPos and not tempNeg:
      signFrict = +1.0
    elif tempNeg and not tempPos:
      signFrict = -1.0
    else:
      print("weird...", tempPos, tempNeg)
      exit()  


    slide = True
    while slide:

        print("  friction direction =", signFrict)
        print("  initial bottom acceleration =", bottomAcceleration(signFrict, finalY))
        if sign(signFrict) != -sign(bottomAcceleration(signFrict, finalY)):   # I check this because I do:  signFrict *= -1.0
            print("aaahhhh!")
            exit()

        # pre calculate
        mu = signFrict * abs(mu_k)
        ratio = (m+M)/M
        ratio2 = mu/ratio
        constTerm = mu * ratio * g_over_r
        constTerm2 = 1 + I/(M*r*r)

        eventTuple = (normalForceSliding, speedDiff)

        tStop = 10/abs(finalY[1])   # feel free to change!
        tList = linspace(finalT, finalT + tStop, num=101)   # nice for making plot; feel free to change num

        solSlide = solve_ivp(derivativesSliding, (finalT, finalT + tStop), (finalY[0],finalY[1],finalY[2],0.0), method = 'Radau', atol = abs_tol, events = eventTuple, t_eval = tList)



        length = len(solSlide.t)
        thetaSlide = solSlide.y[0] % radians(360)
        alphaSlide = [derivativesSliding(solSlide.t[i], solSlide.y[:,i])[1] for i in range(length)]
        FnSlide = [normalForceSliding(solSlide.t[i], solSlide.y[:,i]) for i in range(length)]
        speedDiffArray = [signFrict * speedDiff(solSlide.t[i], solSlide.y[:,i]) for i in range(length)]


        # update graph
        ax[0,1].plot(solSlide.t, thetaSlide, 'k-', solSlide.t, solSlide.y[1], 'r--', solSlide.t, alphaSlide, 'bs', solSlide.t, FnSlide, 'g^', solSlide.t, speedDiffArray, 'm+')
        ax[0,1].set(xlabel="t")
        ax[0,1].legend(['θ mod 2π','ω','α','F_n','r ω + v_x'])


        finalT_hop = False
        if solSlide.t_events[0].size:
            finalT_hop = solSlide.t_events[0][0]
            finalY_hop = solSlide.y_events[0][0]

        finalT_static = False
        if solSlide.t_events[1].size:
            finalT_static = solSlide.t_events[1][0]
            finalY_static = solSlide.y_events[1][0]

        # sort out if you hop, go back to static friction, or neither
        back = False
        if finalT_static and finalT_hop:
            if finalT_static < finalT_hop:   # what if they are equal?
              finalT = finalT_static
              finalY = finalY_static
              back = True
            else:
              finalT = finalT_hop
              finalY = finalY_hop
        elif finalT_static:
            finalT = finalT_static
            finalY = finalY_static
            back = True
        elif finalT_hop:
            finalT = finalT_hop
            finalY = finalY_hop
        else:   # neither
            print("  Make tStop in rollWithSliding() longer.")
            exit()

        if solSlide.status == -1:
          print(" status is -1!")


        print("  final t =", finalT)
        print("  final θ =", str(finalY[0]) + ", so " + str(degrees(finalY[0])%360) + "°")
        print("  final ω =", finalY[1])
        print("  final α =", derivativesSliding(finalT, finalY)[1])
        print("  final F_n =", normalForceSliding(finalT, finalY))
        print("  final v_x =", finalY[2])
        print("  final r ω + v_x =", signFrict * speedDiff(finalT, finalY) )    # if not 0, should have sign opposite signFrict
        print("  mechanical energies:", energyFunc_notHopping(solSlide.y[:,0]), energyFunc_notHopping(finalY))
        print("  work by kinetic friction =", finalY[3])

        slide = False
        if back:
          temp = maxFriction_minus_requiredFriction(finalT, finalY)
          print("  rolling without sliding μ_s F_n - |F_fr| =", temp)
          if temp < 0.0:
            print("  The friction direction changes!")
            slide = True
            signFrict *= -1.0
          else:
            print("  Back to static friction!")
        else:
          print("  Hop!")
        print()



        big[0] = concatenate(( big[0], solSlide.t ))     # t
        big[1] = concatenate(( big[1], solSlide.y[0] ))  # theta
        big[2] = concatenate(( big[2], solSlide.y[1] ))  # omega
        big[3] = concatenate(( big[3], alphaSlide ))     # alpha
        big[4] = concatenate(( big[4], solSlide.y[2] ))  # v_x




    return (finalT, finalY[0:3], back)





#############################
# hop
#############################


def hop(finalT, finalY):

    # for center of mass (CM) initial velocity
    # Note that omega does not change from just before to just after the hop
    vx0 = - r_cm * finalY[1] * sin(finalY[0]) + finalY[2]
    vy0 = r_cm * finalY[1] * cos(finalY[0])



    # y_center is trivially 0 when t=0 and is zero when hoop hits the ground
    # I do not believe that there are ever more solutions than this
    def y_center(t):
      return r_cm * sin(finalY[0]) + vy0 * t - 0.5 * g * t**2 - r_cm * sin(finalY[0] + finalY[1] * t)

    def ay_center(t):
      return -g + r_cm * finalY[1] * finalY[1] * sin(finalY[0] + finalY[1] * t)
    # so that the time derivative of ay_center at t=0 is positive, the hop cannot occur after M points straight up

    def energyFunc_Air(t):
      return 0.5 * (m + M) * (vx0 * vx0 + (vy0 - g * t)**2 ) + weight * (r_cm * sin(finalY[0]) + vy0 * t - 0.5 * g * t**2) + 0.5 * I_cm * finalY[1] * finalY[1]



    temp = ay_center(0)
    print("  initial ay_center =", temp )  # should not be negative
    if not isclose(temp, 0.0, atol = abs_tol):   # I think it even has to be 0!?
      print("oh no!")
      exit()



    # find t, the time the ball is in the air
    bestGuess = 2*vy0/g  # completely ignores rotation of the hoop
    guesses = linspace(3*bestGuess, 0.0, num = 101)   # feel free to change
    absoluteTolerance = 1E-5                          # feel free to change
    for guess in guesses:
      t = fsolve(y_center, guess)[0]
      if not isclose(t, 0.0, atol = absoluteTolerance):
        break
    if isclose(t, 0.0, atol = absoluteTolerance):
      print("  ERROR: fsolve() failed to find non-trivial solution")


    # make graph
    tList = linspace(0.0, 1.25*t, num=101)      # feel free to change num and the final time
    ax[1,1].plot(tList, y_center(tList), 'r-', t, y_center(t), 'r.', tList, 0.0*tList, 'k-')
    ax[1,1].set(xlabel="t_air", ylabel="y_center")


    finalTheta = finalY[0] + finalY[1] * t

    print("  time in air =", str(t) + ", guess =", bestGuess)
    print("  hits ground at t =", finalT + t)
    print("  final θ =", str(finalTheta) + ", so " + str(degrees(finalTheta)%360) + "°")
    print("  final v_x =", vx0 + r_cm * finalY[1] * sin(finalTheta) )
    print("  final r ω + v_x =", vx0 + r_cm * finalY[1] * sin(finalTheta) + r*finalY[1] )
    print("  energies in air:", energyFunc_Air(linspace(0, t, num=3)))
    print()



    tList = linspace(0.0, t, num=101)    # feel free to change num
    omegaAir = tList*0.0 + finalY[1]
    thetaAir = finalY[0] + omegaAir*tList

    big.append( concatenate(( big[0]*0.0, y_center(tList) )) )   # y_pos
    big[0] = concatenate(( big[0], finalT + tList ))
    big[1] = concatenate(( big[1], thetaAir ))
    big[2] = concatenate(( big[2], omegaAir ))
    big[3] = concatenate(( big[3], 0.0*tList ))                  # alpha
    big[4] = concatenate(( big[4], vx0 + r_cm * finalY[1] * sin(finalY[0] + finalY[1] * tList) ))   # v_x

    finalT += t
    # finalY is now [ theta, omega, vx_cm, vy_cm ]
    finalY = [ finalTheta, finalY[1], vx0, vy0 - g*t ]

    return (finalT, finalY) 




#############################
# land
#############################
#
# Assume that the time of the collision (the time until vy_center is 0) is negligibly small.
#
# Energy confirmation is not confirmed here




def trySlide(mu_k):   # +mu_k or -mu_k

    s = sin(finalY[0])
    c = cos(finalY[0])

    term = r_cm*(m+M)/I_cm*(mu_k*(r/r_cm + s) - c)    # a part of the omega calculation
    omega = (finalY[3]*term - finalY[1])/(c*r_cm*term - 1.0)

    vy_cm = r_cm * omega * c
    vx_cm = finalY[2] + mu_k * (vy_cm - finalY[3])
    v_bottom = vx_cm + r_cm * omega * s + r * omega
    v_bottom0 = finalY[2] + r_cm * finalY[1] * s + r * finalY[1]

    # reality check for my sanity; all printed numbers should be 0
    '''
    Fnt = (M + m)*(-finalY[3] + vy_cm)
    Ffrt = (M + m)*(-finalY[2] + vx_cm)
    print( Ffrt - Fnt*mu_k, (r/r_cm + s)*Ffrt - Fnt*c - I_cm/r_cm*(omega-finalY[1]) )
    '''

    # find F_n just after landing
    mu = mu_k
    ratio = (m+M)/M
    ratio2 = mu/ratio
    constTerm = mu * ratio * g_over_r
    constTerm2 = 1 + I/(M*r*r)
    temp = omega*omega*s
    alpha = (-mu*temp + constTerm - g_over_r*c + mu*g_over_r*s - ratio2*temp*s + temp*c/ratio) / (constTerm2 - mu*c - ratio2*s*c - s*s/ratio)
    Fn_after = M_times_r * (alpha * c - temp) + weight

    # (omega, vx_cm, did v_bottom stay same direction?, hop? )
    return ( omega, vx_cm, sign(v_bottom)==sign(v_bottom0), Fn_after < 0.0 )



def tryStopBottom():

    # pre calculate
    c = cos(finalY[0])
    s = sin(finalY[0])
    term1 = r/r_cm + s
    term2 = r_cm * (m+M) / I_cm

    # get omega
    omega = (finalY[1] + c * term2 * finalY[3] - term1 * term2 * finalY[2]) / (1.0 + c**2 * term2 * r_cm + term1**2 * term2 * r_cm)

    vx_cm = -r_cm * omega * s - omega * r

    # reality check for my sanity
    '''
    vy_cm = r_cm * omega * c
    Fnt = (M + m)*(-finalY[3] + vy_cm)
    Ffrt = (M + m)*(-finalY[2] + vx_cm)
    print( term1*Ffrt - Fnt*c - I_cm/r_cm*(omega-finalY[1]) )    # should be 0
    '''

    # find F_n just after landing
    alpha = - (omega*omega + g_over_r) * c / (aTerm + 2*s)
    Fn_after = M_times_r * ( alpha * c - omega*omega*s ) + weight

    # find static F_fr after landing
    Ffr_after = - M_times_r * ( omega*omega*c + alpha*(1+s) ) - m_times_r * alpha

    # (omega, vx_cm, static after landing?, hop? )
    return ( omega, vx_cm, mu_s * Fn_after > abs( Ffr_after ), Fn_after < 0.0 )



def trySlideAgain(mu_k):   # +mu_k or -mu_k is the initial friction direction but then friction changes direction

    # pre calculate
    s = sin(finalY[0])
    c = cos(finalY[0])
    term1 = r/r_cm + s
    term2 = r_cm * (m+M) / I_cm
    ratioA = ( mu_k + term1*term2*r_cm*(mu_k*term1-c) ) / ( -mu_k + term1*term2*r_cm*(-mu_k*term1-c) )

    v_bottom0 = finalY[2] + r_cm * finalY[1] * s + r * finalY[1]

    # get frac, the fraction of Fn's impulse that gives the original friction direction
    # A, B, and C are for quadratic equation
    A = 2*mu_k*ratioA*( c*finalY[1]*r_cm - finalY[3] + term1*term2*(c*r_cm*v_bottom0/ratioA - c*finalY[2]*r_cm - finalY[3]*r - finalY[3]*r_cm*s) )
    B = mu_k*term1*term2 * (c*finalY[2]*r_cm - 3*c*r_cm*v_bottom0/ratioA + finalY[3]*r + finalY[3]*r_cm*s ) + c*r_cm*term2 * (c*finalY[2] - c*v_bottom0/ratioA + finalY[3]*r/r_cm + finalY[3]*s) - c*finalY[1]*mu_k*r_cm + finalY[1]*r + finalY[1]*r_cm*s + finalY[2] + finalY[3]*mu_k - v_bottom0/ratioA
    B *= ratioA
    C = v_bottom0 * ( c*r_cm*term2 * ( c + mu_k*term1 ) + 1.0)
    frac1 = (-B + sqrt(B**2 - 4*A*C))/(2*A)
    frac2 = (-B - sqrt(B**2 - 4*A*C))/(2*A)
    if 0 < frac1 < 1:
      if 0 < frac2 < 1:
        return (False,False,False,False)  # bazoinga!! Both work!
      frac = frac1
    elif 0 < frac2 < 1:
      frac = frac2
    else:
      return (False,False,False,False)  # bazoinga!! Neither work!

    # get omega
    term = term2*(mu_k*term1*(2*frac-1) - c)    # a part of the omega calculation
    omega = (-finalY[1] + finalY[3]*term)/(-1.0 + c*r_cm*term)

    vy_cm = r_cm * omega * c
    vx_cm = finalY[2] + mu_k * (2*frac - 1) * (vy_cm - finalY[3])
    v_bottom = vx_cm + r_cm * omega * s + r * omega

    # reality check for my sanity; all printed numbers should be 0
    '''
    Fnt = (M + m)*(-finalY[3] + vy_cm)   # an impulse
    Ffrt = (M + m)*(-finalY[2] + vx_cm)  # an impulse
    print( Ffrt - Fnt*mu_k*(2*frac - 1), Ffrt*term1 - Fnt*c - I_cm*(-finalY[1] + omega)/r_cm, frac*ratioA*v_bottom + (1 - frac)*v_bottom0 )
    '''

    # find F_n just after landing
    mu = -mu_k
    ratio = (m+M)/M
    ratio2 = mu/ratio
    constTerm = mu * ratio * g_over_r
    constTerm2 = 1 + I/(M*r*r)
    temp = omega*omega*s
    alpha = (-mu*temp + constTerm - g_over_r*c + mu*g_over_r*s - ratio2*temp*s + temp*c/ratio) / (constTerm2 - mu*c - ratio2*s*c - s*s/ratio)
    Fn_after = M_times_r * (alpha * c - temp) + weight

    # (omega, vx_cm, did v_bottom change directions?, hop? )
    return ( omega, vx_cm, sign(v_bottom) == -sign(v_bottom0), Fn_after < 0.0 )





def land(finalY):

    v_bottom0 = finalY[2] + r_cm * finalY[1] * sin(finalY[0]) + r * finalY[1]
    signFrict = -sign(v_bottom0)   # initially at least

    trySlideSimple = trySlide(signFrict * abs(mu_k))
    tryStop = tryStopBottom()
    trySlideComplex = trySlideAgain(signFrict * abs(mu_k))

    static = False
    print("  Landing!")
    if trySlideSimple[2]:
      print("  kinetic friction direction =", signFrict)
      omega, vx_cm, _, hop = trySlideSimple
    elif tryStop[2]:
      print("  v_bottom becomes 0")
      omega, vx_cm, _, hop = tryStop
      static = True
    elif trySlideComplex[2]:
      print("  kinetic friction initial but not final direction =", signFrict)
      omega, vx_cm, _, hop = trySlideComplex
    else:
      print("ERROR")
      print(trySlideSimple)
      print(tryStop)
      print(trySlideComplex)
      exit()


    finalY = [finalY[0], omega, vx_cm + r_cm*omega*sin(finalY[0])]

    if hop:
      print("  another hop occurs before vy_center can become zero!?")
      print()
      return (finalY, static, hop)


    print("  final ω =", omega)
    print("  final v_x =", finalY[2])

    # Should be 0 if and only if v_bottom = 0 solution was used
    # If not 0, should have sign opposite sign of friction
    print("  final r ω + v_x =", r*omega + finalY[2] )

    print("  energy =", energyFunc_notHopping(finalY))
    if static:
      print("  Rolling without sliding!")
    else:
      print("  Sliding!")
    print()

    return (finalY, static, hop)



#############################
# driver code
#############################



# Require rolling without sliding from initial conditions
print()
if normalForce(0, (theta0,omega0)) <= 0.0:
  print("  Hopping at t=0!")
  exit()
if maxFriction_minus_requiredFriction(0, (theta0,omega0)) <= 0.0:
  print("  Sliding at t=0!")
  exit()



# big is a 2D list that continuously grows as the code executes
#   [ time, theta, omega, alpha, v_x ]
big = [ [] for i in range(5) ]


# initialize
finalT = 0.0
finalY = (theta0, omega0)

while True:

  finalT, finalY, slide = rollWithoutSliding(finalT, finalY[0:2])   # will exit the code if nothing happens before tStop

  if slide:
    finalT, finalY, back = rollWithSliding(finalT, finalY)   # will handle sign of friction changing directions
  else:      # hop
    break

  if not back:   # hop
    break


finalT, finalY = hop(finalT, finalY)

land(finalY)


makePlots()
animate()


