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
            //Sample angle from 0 to 2PI
            int numPoints = 1001;
            double[] phiArray = new double[numPoints];
            double step = Math.PI / (numPoints - 1); // subtract 1 to include both endpoints

            for (int i = 0; i < numPoints; i++)
            {
                phiArray[i] = i * step;
            }

            //Loop to find the max value of Klein Nichina Cross Section
            foreach (double phi in phiArray)
            {
                
                //Calculate energy of the scattere photon
                var photonEneryScat = photonEnery / (1 + (photonEnery / 0.511) * (1 - Math.Cos(phi)));
                var thisKlNichCrs = (Math.PI) * (Math.Sin(phi)) * (Math.Pow((photonEneryScat / photonEnery), 2));
                thisKlNichCrs = thisKlNichCrs * (photonEnery / photonEneryScat + photonEneryScat / photonEnery - Math.Pow(Math.Sin(phi),2));
                thisKlNichCrs = thisKlNichCrs * Math.Pow(electronRadius,2);

                if (thisKlNichCrs > maxKlNich) 
                    maxKlNich = thisKlNichCrs;
            }
            return maxKlNich;
        }
    }
}
