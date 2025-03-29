using MCassignment_PHYS539_PeidusD_2025.Helpers;

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public static class maxCollisionScatteringProbability
    {
        public static double GetMaxColScatteringProb(double eMeV,double pathLength)
        {
            var thetaSquared = PhysicsData.GetMeanSquaredScatteringAngleSquaredForWater(eMeV, pathLength);
            var maxHardCollisionScatteringProbabilityValue = 0.0;
            var numberOfSteps = 1001;
            var maxAngle = 0.0;
            double[] stepsSampling = new double[numberOfSteps];
            var step = 0.001;
            for (int i = 0; i < numberOfSteps; i++)
            {
                stepsSampling[i] = i * step;
            }

            foreach(var angleRad in stepsSampling)
            {
                var a = 2*angleRad/Math.Pow(Math.PI,2); 
                a = a / thetaSquared;

                var b = Math.Pow(angleRad, 2); 
                b = b / thetaSquared; 
                b = Math.Exp(-b);

                var thisScatteringProbability = a * b;

                if (thisScatteringProbability > maxHardCollisionScatteringProbabilityValue)
                { maxHardCollisionScatteringProbabilityValue = thisScatteringProbability;
                    maxAngle = angleRad;
                }
            }
            return maxHardCollisionScatteringProbabilityValue;
        }
        
    }
}
