

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public static class PairProductionHelper
    {
        //Determine max possible Klein-Nichina cross-section using the rejection method
        public static double SamplePairProductionCrSct(double particleEnergy)
        {
            var maxPairCs = 0.0;
            //Sample angle from 0 to PI
            int numPoints = 101;
            double[] phiArray = new double[numPoints];
            double step = Math.PI / (numPoints - 1); // subtract 1 to include both endpoints

            for (int i = 1; i < numPoints; i++)
            {
                phiArray[i] = i * step;
            }
            foreach (double phi in phiArray)
            {
                double beta = Math.Sqrt(1 - 1 / Math.Pow((particleEnergy / 0.511) + 1, 2));

                var thisPairProdCs = (Math.PI * Math.Sin(phi)) / Math.Pow((1 - beta * Math.Cos(phi)), 2);

                if (thisPairProdCs > maxPairCs)
                {
                    maxPairCs = thisPairProdCs;
                }
            }
            return maxPairCs;
        }
    }
}
