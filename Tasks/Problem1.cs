using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Engine;
using MCassignment_PHYS539_PeidusD_2025.Helpers;
using ScottPlot;
namespace MCassignment_PHYS539_PeidusD_2025.Tasks
{
    public static class Problem1
    {
        /************************************************************************************
        Class is a controller fot the Problem I. Monte Carlo simulation of a pencil beam.
        Photons of 3 energies are simulated. Output is designed specifically to fulfill
        the requirements of the Assignment #2 PHYS539. I did not develop it for any other purpose,
        so it outputs .CSV files with the peak Y slice for the projection XZ for isolines plotting. 
        Isolines are then plotted with using an external software. Also, for functionality 
        assessment, statistics is collected into a separate container, which allows manual data
        review and export.
        ************************************************************************************/
        public static void KernelGenerationTask(StatisticsToExport statistics)
        {
            //Array of energies of interest.
            double[] photonEnergy = [2.0, 6.0, 10.0];

            //Simulate 10e6 particles.
            int numPhotons = 1000000;

            //Create 12.5cm3 phantom.
            double voxelSizeCm = 0.05;
            var nx = 250;
            var ny = 250;
            var nz = 250;

            //Main loop foe arch energy
            foreach (var _photonEnergy in photonEnergy)
            {
                Console.WriteLine($"\n{_photonEnergy}MeV photon simulation started...");

                //Sample max KN probability distribution value for furthe use in the Rejection Method.
                var maxKlNichProb = maxKleinNichinaHelper.SampleMaxKleinNichina(_photonEnergy);

                //Create a container of energy values within the phantom.
                var voxels = new WaterTankPhantom_UniformResolution(nx, ny, nz, voxelSizeCm);

                //Sample the lookup table for defining scattering angles for charged particles.
                Dictionary<double, double> lookupElScat = new();
                lookupElScat= ElectronPositronScatteringAngleLookupTable.LookupElectronScatAngle(_photonEnergy, voxelSizeCm);

                //Call a method to initiate the simulation.
                var engine = new MonteCarloEngine();
                engine.RunPhotonKernelCalculation(voxels, numPhotons, _photonEnergy, maxKlNichProb, statistics,voxelSizeCm, lookupElScat);

                //Define .CSV file name
                string baseName = $"Kernel_{(double)_photonEnergy}MeV.CSV";

                //voxels.ExportFullKernel(baseName, numPhotons, _photonEnergy);
                Console.WriteLine($"Kernel generation complete for {(double)_photonEnergy}MeV");

                //Metrics to check functionality
                var data = statistics.GetPhotoelectricScatteringAngles();
                double meanComptonAngle = statistics.AverageAngleCompton();
                double meanPairAngle_e = statistics.AverageAnglePair_e();
                double meanPairAngle_p = statistics.AverageAnglePair_p();
                double meanHardCollisionAngle = statistics.AverageHardCollisionAngle();

                // Export statistics
                var filePath = "Statistics.csv"; // Add .csv extension
                using (var writer = new StreamWriter(filePath))
                {
                    // Optional: Write header
                    writer.WriteLine("Statistics");

                    // Write each value in a new line
                    int iterator = 0;
                    foreach (var angle in StatisticsToExport.scatteringAnglesComptonPhotons)
                    {
                        iterator = iterator + 1;
                        writer.WriteLine(angle.ToString());
                        if (iterator > 200000) break; 
                        //We do not neeed all 10e6 points for reporting. 
                        //Limiting this number helps to avoid Excel crashes. 
                        //My laptop needs upgrade...
                    }
                }

                Console.WriteLine("CSV export completed: " + filePath);

                //Clear statistics before commencing on the new energy.
                statistics.ClearComptonEnergy();
                statistics.ClearCompton();
                statistics.ClearPair_e();
                statistics.ClearPair_p();
                statistics.ClearHardCollisionAngle();
            }
        }
    }
}
