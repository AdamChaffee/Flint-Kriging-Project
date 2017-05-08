# Flint-Kriging-Project

## Overview and History
The Flint water crisis is an environmental catastrophe that has garnered media attention over the last few years. In 2014 local officials changed the City of Flint, Michigan's water source from Lake Huron to the nearby Flint river in an attempt to save money. After making the switch, government officials failed to apply corrosion inhibitors to the water before piping it into the city system.

Much of Flint's water infrastructure was built in the early 1900's at a time when lead and copper lining of water pipes was common. The Flint river water irreversibly corroded the pipes causing lead and copper linings to leach into the water supply. Thousands of Flint residents who drank the water suffered a range of health problems. Children were especially vulnerable as they became susceptible to physical growth, behavior, and learning problems.

After switching water sources it quickly became apparent that water quality had diminished. The water became brown and murky, and several residents and their children began to experience significant health issues. State officials came in to test the water in a variety of locations and found no cause for mandatory action.

Local residents and scientists were not convinced. A grant from the MacArthur Foundation put testing kits and instructions into the hands of local residents who faithfully executed testing and reporting instructions to the best of their ability.

## The data
The data I used comes from resident samples of water from their own homes from January to March, 2017. Up until recently this data has been controversial as it does not come from scientific experts. However, recent research has demonstrated that citizen testing is more accuracte and informative than the state testing due to more comprehensive sampling and better adherence to the testing instructions. For more information, the April 2017 edition of "Significance", a monthly statistics publication, contains a great story about this controversy.

## Research goals
The best way to identify potential lead leaching is to look at city construction records. Unfortunately many of the records are lost and over half of the city contains water pipes of unknown composition. I realized statistics could be used to help fill in the gaps.

I decided to utilize several kriging methods using residential testing data in an attempt to identify potential regions of high lead leaching. I explored ordinary, simple, universal, and indicator kriging to obtain estimates of lead leaching levels throughout the city. Since the testing data also contained copper levels, I explored using copper as a covariate in co-kriging as well.

## Results
After using leave-one-out cross validation on all kriging methods, co-kriging with copper produced by far the lowest residual sum of squares estimate due to the weak positive relationship between copper and lead. When comparing the contour map to construction records, the areas of concern my map identified line up well with known areas constructed with lead pipe linings. However, no new areas of concern were identified from these methods.

I presented this work in March, 2017. Presentation slides are available upon request.

## References
Data obtained from: http://www.michigan.gov/flintwater/0,6092,7-345-76292_76294_76297---,00.html

Article from Significance Magazine: http://onlinelibrary.wiley.com/doi/10.1111/j.1740-9713.2017.01016.x/full
