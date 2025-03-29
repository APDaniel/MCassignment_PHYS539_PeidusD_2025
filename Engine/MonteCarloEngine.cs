using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Helpers;
using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Engine
{
    public class MonteCarloEngine
    {
        public void RunPhotonKernelCalculation(WaterTankPhantom_UniformResolution voxels, 
            int numPhotons, 
            double photonEnergy, 
            double maxKleinNichProb, 
            StatisticsToExport statistics, double pathlength, 
            Dictionary<double, double> lookupElScat)

        {
            for (int i =1; i <= numPhotons; i++)
            {
                Random random = new Random();
                //Launch photons from center of XY at Z=0
                double centerX = (voxels.Nx * voxels.VoxelSizeCm) / 2;
                double centerY = (voxels.Ny * voxels.VoxelSizeCm) / 2;
                double centerZ = (voxels.Nz * voxels.VoxelSizeCm) / 10;
                Vector3 position = new Vector3((float)centerX, (float)centerY, (float)centerZ);
                Vector3 direction = new Vector3(0f, 0f, (float)voxels.VoxelSizeCm);
                SimulatePhotonInteraction(voxels, position,direction, photonEnergy, maxKleinNichProb, statistics, pathlength, lookupElScat);
            }

            int peakY = voxels.FindPeakYSlice();
            voxels.SaveSliceWithTracksAsImage(peakY, "kernel_peak_with_tracks.png", allTracks, photonEnergy);
            voxels.ExportSliceY(peakY,"kernexMaxSlice.csv",numPhotons,photonEnergy);
        }

        public void RunPhotonKernelCalculationWithTravelDistance(WaterTankPhantom_UniformResolution voxels,
            int numPhotons,
            double photonEnergy,
            double maxKleinNichProb,
            StatisticsToExport statistics, double pathlength,
            Dictionary<double, double> lookupElScat)

        {
            var phPerVxl = numPhotons/(1600);
            for (double ix = 0; ix < voxels.Nx; ix++)
            {
                if (ix >= 75 && ix <= 125)
                {
                    for (double iy = 0; iy < voxels.Ny; iy++)
                    {
                        if (iy >= 75 && iy <= 125)
                        for (int photon = 1; photon <= phPerVxl; photon++)
                        {
                            Random random = new Random();
                            var R = random.NextDouble();
                            var iz = -1 / 0.022;
                            R = Math.Log(R);
                            iz = iz * R;
                            iz = iz / voxels.VoxelSizeCm;
                            //Launch photons for each voxel equally on surface and modeulate depth of penetration
                            Vector3 position = new Vector3(
                                (float)ix* (float)voxels.VoxelSizeCm, 
                                (float)iy*(float)voxels.VoxelSizeCm, 
                                (float)iz);

                            Vector3 direction = new Vector3(0f, 0f, (float)voxels.VoxelSizeCm);
                            SimulatePhotonInteraction(voxels, position, direction, photonEnergy, maxKleinNichProb, statistics, pathlength, lookupElScat);
                        }
                    }
                }  
            }


            int peakY = voxels.FindPeakYSlice();
            voxels.SaveSliceWithTracksAsImage(peakY, "kernel_peak_with_tracks.png", allTracks, photonEnergy);
            voxels.ExportSliceY(peakY, "kernexMaxSlice.csv", numPhotons, photonEnergy);
        }
        public void RunElectronBroadBeamSimulation(WaterTankPhantom_ReducedResolution voxels, 
            double electronsEnteredPerVoxel, double electronEnergy, 
            StatisticsToExport statistics, double numElectrons)
        {
            MonteCarloEngine engine = new();
            var nx = voxels.Nx;
            var ny = voxels.Ny;
            var nz = voxels.Nz;

            for (double i = 0; i < ny - 1; i=i+1)
            {
                for (double i1 = 0; i1 < nx - 1; i1=i1+1)
                {
                    for (int i2 = 0; i2 < electronsEnteredPerVoxel; i2++)
                    {
                        var coordinateX = i1 * voxels.VoxelSizeCm_XY; 
                        var coordinateY = i * voxels.VoxelSizeCm_XY;
                        var stepSize_Z = voxels.VoxelSizeCm_Z;

                        if (coordinateX>=4&& coordinateX <= 8)
                        {
                            if (coordinateY >= 4 && coordinateY <= 8)
                            {
                                Vector3 startPos = new((float)coordinateX, (float)coordinateY, 0);
                                Vector3 startDirection = new(0f, 0f, (float)voxels.VoxelSizeCm_Z);
                                Vector3 deflectionOnTheSurface = new();
                                var theta=PhysicsData.SampleMultipleScatteringAngleTranslation(electronEnergy, voxels.VoxelSizeCm_Z);
                                var phi = PhysicsData.SampleAzimutalAngle();
                                Vector3 scatterDirection = PhysicsData.SphericalToCartesian(theta, phi);

                                //Switch to on/off mult scattering (also, need to turn on/of in the particle transport method
                                //deflectionOnTheSurface = VectorRotationHelper.RotateToDirection(scatterDirection, startDirection);
                                deflectionOnTheSurface = startDirection;

                                Dictionary<double, double> empty = new();
                                TransportChargedParticles_ReducedRes(voxels, startPos,
                                    electronEnergy, deflectionOnTheSurface,
                                    electronEnergy, deflectionOnTheSurface,
                                    0.1, statistics, empty, voxels.VoxelSizeCm_Z);
                            }
                        }
                    }
                }
            }

            int peakY = voxels.FindPeakYSlice();
            voxels.SaveSliceWithTracksAsImage(peakY, "kernel_peak_with_tracks.png", allTracks, 10);
            voxels.ExportSliceY(peakY, "10MeV_beta-_kernexMaxSlice.csv", numElectrons, 10);

        }
        private static void SimulatePhotonInteraction(
            WaterTankPhantom_UniformResolution voxels, 
            Vector3 position, Vector3 direction,
            double energy, double maxKleinNichProb,
            StatisticsToExport statistics, double electronPositronPathLength,
            Dictionary<double, double> lookupElScat)
        {
            //Sample interaction point and direction
            Vector3 interactionPoint = position + direction;

            //Select interaction type
            string interaction = PhysicsData.SampleInteractionType(energy);

            //Compute electron parameters and transport it
            var chargedParticleFromInteraction = PhysicsData.GetElectronFromInteraction(interaction,energy,direction,statistics, maxKleinNichProb);
            chargedParticleFromInteraction.positronDirection.X = -chargedParticleFromInteraction.positronDirection.X;
            chargedParticleFromInteraction.positronDirection.Y = -chargedParticleFromInteraction.positronDirection.Y;

            //Transport charged particles until their energy is depleted.
            //e- deposits its rest energy when stopped.
            //e+ do not due to annihilation.
            TransportChargedParticles(voxels,interactionPoint, 
                chargedParticleFromInteraction.electronEnergy, chargedParticleFromInteraction.electronDirection,
                chargedParticleFromInteraction.positronEnergy, chargedParticleFromInteraction.positronDirection,
                electronPositronPathLength,statistics, lookupElScat);
        }
        private static List<List<Vector3>> allTracks = new();
        public static void TransportChargedParticles(
            WaterTankPhantom_UniformResolution voxels, Vector3 startPos, 
             double electronEnergy, Vector3 electronDirection,
             double positronEnergy, Vector3 positronDirection,
             double pathLength, StatisticsToExport statistic,
             Dictionary<double, double> lookupElScat)
        {
            const double stepSize = 0.05; //in cm
            const double minEnergy = 0.02; //MeV cutoff

            Vector3 position = startPos;
            var currentTrack = new List<Vector3>();
            currentTrack.Add(position);

            //Transport electron and deposit its rest energy when it stopped
            if (electronEnergy > 0)
            {
                while (electronEnergy > minEnergy && voxels.BoundaryCheck(position.X, position.Y, position.Z))
                {
                    //Get stopping power from the CSDA range
                    double stoppingPower = CSDArangeHelper.GetStoppingPower(electronEnergy);

                    //Simulate energy loss for the particle
                    double energyLoss = Math.Min(electronEnergy, stoppingPower * stepSize);

                    //Simulate multiple scattering accounting for elastic scattering
                    //double theta = PhysicsData.SampleMultipleScatteringAngle(electronEnergy, pathLength);

                    double theta = 0.0;
                    /*
                    foreach (var key in lookupElScat.Keys)
                    {
                        if (key > electronEnergy)
                        {
                            theta = lookupElScat[key];
                            break;
                        }      
                    }
                    */
                    theta = PhysicsData.SampleMultipleScatteringAngleTranslation(electronEnergy, pathLength);
                    electronEnergy -= energyLoss;

                    statistic.RecordHardCollisionAngle(theta);
                    double phi = PhysicsData.SampleAzimutalAngle();
                    Vector3 scatterDirection = PhysicsData.SphericalToCartesian(theta, phi);
                    electronDirection = VectorRotationHelper.RotateToDirection(scatterDirection, electronDirection);

                    //Take a step
                    position += electronDirection * (float)stepSize;
                    currentTrack.Add(position);

                    //Deposit energy into voxel of interaction
                    voxels.DepositEnergy(position.X, position.Y, position.Z, energyLoss);
                }
                //Deposit rest energy into the very last voxel of e- track
                voxels.DepositEnergy(position.X, position.Y, position.Z, 0.511);
            }

            //Check if the electron was produced with a positron of non-zero energy and transport it
            if (positronEnergy > 0)
            {
                while (positronEnergy > minEnergy && voxels.BoundaryCheck(position.X, position.Y, position.Z))
                {
                    //Get stopping power from the CSDA range
                    double stoppingPower = CSDArangeHelper.GetStoppingPower(positronEnergy);

                    //Simulate energy loss for the particle
                    double energyLoss = Math.Min(positronEnergy, stoppingPower * stepSize);
                    positronEnergy -= energyLoss;

                    //Simulate multiple scattering
                    //double theta = PhysicsData.SampleMultipleScatteringAngle(electronEnergy, pathLength);
                    double theta = PhysicsData.SampleMultipleScatteringAngleTranslation(positronEnergy,pathLength);
                    /*
                    foreach (var key in lookupElScat.Keys)
                    {
                        if (key > electronEnergy)
                        {
                            theta = lookupElScat[key];
                            break;
                        }
                    }
                    */
                    double phi = PhysicsData.SampleAzimutalAngle();
                    Vector3 scatterDirection = PhysicsData.SphericalToCartesian(theta, phi);
                    positronDirection = VectorRotationHelper.RotateToDirection(scatterDirection, positronDirection);

                    //Take a step
                    position += positronDirection * (float)stepSize;
                    currentTrack.Add(position);

                    //Deposit energy into voxel of interaction
                    voxels.DepositEnergy(position.X, position.Y, position.Z, energyLoss);
                }
            }
            allTracks.Add(currentTrack);

        }

        public static void TransportChargedParticles_ReducedRes(
            WaterTankPhantom_ReducedResolution voxels, Vector3 startPos,
             double electronEnergy, Vector3 electronDirection,
             double positronEnergy, Vector3 positronDirection,
             double pathLength, StatisticsToExport statistic,
             Dictionary<double, double> lookupElScat, double stepSize)
        {
            const double minEnergy = 0.02; //MeV cutoff

            Vector3 position = startPos;
            var currentTrack = new List<Vector3>();
            currentTrack.Add(position);

            //Transport electron and deposit its rest energy when it stopped
            if (electronEnergy > 0)
            {
                while (electronEnergy > minEnergy && voxels.BoundaryCheck(position.X, position.Y, position.Z))
                {
                    //Get stopping power from the CSDA range
                    double stoppingPower = CSDArangeHelper.GetStoppingPower(electronEnergy);

                    //Simulate energy loss for the particle
                    double energyLoss = Math.Min(electronEnergy, stoppingPower * stepSize);

                    //Simulate multiple scattering accounting for elastic scattering
                    double theta=0.0;
                    /*
                    foreach (var key in lookupElScat.Keys)
                    {
                        if (key > electronEnergy)
                        {
                            theta = lookupElScat[key];
                            break;
                        }      
                    }
                    */
                    //theta = PhysicsData.SampleMultipleScatteringAngleTranslation(electronEnergy, pathLength);
                    electronEnergy -= energyLoss;

                    statistic.RecordHardCollisionAngle(theta);
                    double phi = 0.0;
                    //double phi = PhysicsData.SampleAzimutalAngle();


                    Vector3 scatterDirection = PhysicsData.SphericalToCartesian(theta, phi);

                    //Can turn on/off for account for multiple scattering
                    electronDirection = VectorRotationHelper.RotateToDirection(scatterDirection, electronDirection); 

                    //Take a step
                    position += electronDirection * (float)stepSize;
                    currentTrack.Add(position);

                    //Deposit energy into voxel of interaction
                    voxels.DepositEnergy(position.X, position.Y, position.Z, energyLoss);
                }
                //Deposit rest energy into the very last voxel of e- track
                voxels.DepositEnergy(position.X, position.Y, position.Z, 0.511);
            }
            allTracks.Add(currentTrack);

        }
    }
}
