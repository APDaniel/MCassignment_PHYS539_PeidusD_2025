using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public static class maxKleinNichinaHelper
    {
        //Determine max possible Klein-Nichina cross-section using the rejection method
        public static double SampleMaxKleinNichina(double photonEnery)
        {
            var maxKlNich = 0.0;
            var electronRadius = Constants.electronRadius_m;
            //Sample angle from 0 to PI
            int numPoints = 1001;
            double[] phiArray = new double[numPoints];
            double step = Math.PI / (numPoints - 1); // subtract 1 to include both endpoints

            for (int i = 1; i < numPoints; i++)
            {
                phiArray[i] = i * step;
            }

            //Loop to find the max value of Klein Nichina Cross Section
            var maxPhy = 0.0;
            foreach (double phi in phiArray)
            {
                //Calculate energy of the scattered photon
                var hv = photonEnery;
                var hv_prime = hv/(1+(hv/0.511)*(1-Math.Cos(phi)));
                var a = Math.PI * Math.Sin(phi) * Math.Pow(electronRadius, 2);
                var b = Math.Pow((hv_prime / hv), 2);
                var c = hv / hv_prime + hv_prime / hv - Math.Pow(Math.Sin(phi), 2);
                var cs = a * b * c;

                if (cs > maxKlNich)
                { 
                    maxKlNich = cs;
                    maxPhy = phi;
                }
            }
            return maxKlNich;
        }
    }
}
