using MCassignment_PHYS539_PeidusD_2025.Helpers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Containers
{
    public class ElectronPositronScatteringAngleLookupTable
    {
        private Dictionary<double, double> _lookUpTable { get; set; }

        public ElectronPositronScatteringAngleLookupTable(Dictionary<double, double> lookupTable)
        {
            _lookUpTable = lookupTable;
        }

        public static Dictionary<double, double> LookupElectronScatAngle(double eMeV, double pathLength) 
        {
            Dictionary<double, double> lu = new();

            var numberOfSamples = eMeV * 100;
            var step = 0.0;
            step = eMeV / numberOfSamples;

            double sampledEnergy=0.0; 
            sampledEnergy = sampledEnergy - step;

            for (int i=0;i<= numberOfSamples; i++)
            {
                sampledEnergy = sampledEnergy + step;
                var scatteringAngle = PhysicsData.SampleMultipleScatteringAngleTranslation(sampledEnergy, pathLength); 
                //PhysicsData.SampleMultipleScatteringAngle(sampledEnergy, pathLength);

                //scatteringAngle = 0.02;
                lu.Add(sampledEnergy, scatteringAngle);
            }
            //Console.WriteLine($"Electron scattering angle energy lookup table calculated for {eMeV} MeV.");
            return lu;
        }

    }
}
