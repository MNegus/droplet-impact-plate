# plate-impact
DNS code for droplet impact onto moving plates using [Basilisk](<http://basilisk.fr/>). 

Contributors: [Radu Cimpeanu](<https://warwick.ac.uk/fac/sci/maths/people/staff/cimpeanu/>), [Matthew Moore](<https://www.maths.ox.ac.uk/people/matthew.moore>), [James Oliver](<https://www.maths.ox.ac.uk/people/james.oliver>), [Ben Fudge & Alfonso Castrejón-Pita](<http://www2.eng.ox.ac.uk/fluidlab/members>)]

## Preamble
* [Basilisk](<http://basilisk.fr/>) is an open source volume of fluid (VOF) package that can be used to solve the incompressible Navier-Stokes equations
* In particular, Basilisk is well suited to solving two-phase flow problems such as droplet impact

## The problem
* We wish to study the dynamics of droplets impacting onto deformable substrates, such as elastic sheets or cling film (see https://aip.scitation.org/doi/full/10.1063/1.2969755).
* As a first step towards this goal, we would like to conduct simulations of impact onto flat plates which are restricted to moving in the vertical direction
* This requires embedding a solid into the computational domain of the problem
* Basilisk currently only supports embedding solids in two-phase problems when only one of the phases makes contact with the solid (such as flow over a bump http://basilisk.fr/src/examples/gaussian-ns.c)
* The issue here is that both the air and the liquid phase will meet the substrate during droplet impact, so we cannot embed a solid in the same way that is currently implemented in Basilisk

## Proposed solution
* We can create artifical solids by using "Stephane's Trick" to move a cylinder through a liquid http://basilisk.fr/sandbox/popinet/movingcylinder.c
* The idea is to create another volume fraction (here we call it plate[]), which is 1 inside the plate and 0 outside of the plate
* In order for the plate to act as a solid, we need to impose the kinematic condition on the velocity (that is the vertical velocity of the fluid needs to match the plate speed)

## Goal
Conduct simulations of axisymmetric droplet impact onto a moving plate, where the displacement of the plate is determined by the forces felt on it from the droplet. In particular, we wish to model the problem of a plate being suspended by a spring and a dashpot, in order to study how elastic and dissipative effects on the movement affect the dynamics of the droplet. 

## Structure of repository
* This repository is to document the progress we are making into studying impact onto moving plates. Each sub directory documents a different task to complete in order to solve the eventual problem of impact onto a plate-spring-dashpot system