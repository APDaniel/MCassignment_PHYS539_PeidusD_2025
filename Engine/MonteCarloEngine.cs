using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Helpers;
using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Engine
{
    public class MonteCarloEngine
    {
        public void RunKernelCalculation(WaterTankPhantom voxels, 
            int numPhotons, 
            double photonEnergy, 
            double maxKleinNichProb, 
            StatisticsToExport statistics)
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
                SimulatePhotonInteraction(voxels, position,direction, photonEnergy, maxKleinNichProb, statistics);
            }
            int peakY = voxels.FindPeakYSlice();
            voxels.SaveSliceWithTracksAsImage(peakY, "kernel_peak_with_tracks.png", allTracks, photonEnergy);

        }
        private static void SimulatePhotonInteraction(
            WaterTankPhantom voxels, 
            Vector3 position, Vector3 direction,
            double energy, double maxKleinNichProb,
            StatisticsToExport statistics)
        {
            //Sample interaction point and direction
            Vector3 interactionPoint = position + direction;

            //Select interaction type
            string interaction = PhysicsData.SampleInteractionType(energy);

            //Compute electron parameters and transport it
            var electronFromInteraction = PhysicsData.GetElectronFromInteraction(interaction,energy,direction,statistics, maxKleinNichProb);

            TransportChargedParticles(voxels,interactionPoint, 
                electronFromInteraction.electronEnergy, electronFromInteraction.electronDirection,
                electronFromInteraction.positronEnergy, electronFromInteraction.positronDirection);
        }
        private static List<List<Vector3>> allTracks = new();
        private static void TransportChargedParticles(
            WaterTankPhantom voxels, Vector3 startPos, 
             double electronEnergy, Vector3 electronDirection,
             double positronEnergy, Vector3 positronDirection)
        {
            const double stepSize = 0.05; //in cm
            const double minEnergy = 0.02; //MeV cutoff

            Vector3 position = startPos;
            var currentTrack = new List<Vector3>();
            currentTrack.Add(position);

            //Transport electron and deposit its rest energy when it stopped
            while (electronEnergy > minEnergy && voxels.BoundaryCheck(position.X, position.Y, position.Z))
            {
                //Get stopping power from the CSDA range
                double stoppingPower = CSDArangeHelper.GetStoppingPower(electronEnergy);

                //Simulate energy loss for the particle
                double energyLoss = Math.Min(electronEnergy, stoppingPower * stepSize);
                electronEnergy -= energyLoss;

                //Simulate multiple scattering accounting for elastic scattering
                double theta = PhysicsData.SampleMultipleScatteringAngle(electronEnergy);
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
                    double theta = PhysicsData.SampleMultipleScatteringAngle(positronEnergy);
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
    }
}
