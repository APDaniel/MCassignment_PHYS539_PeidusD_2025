using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Helpers;
using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Engine
{
    public class MonteCarloEngine
    {
        public void Run(WaterTankPhantom voxels, int numPhotons, double photonEnergy)
        {
            for (int i =1; i <= numPhotons; i++)
            {
                SimulatePhoton(voxels, photonEnergy);
            }
            int peakY = voxels.FindPeakYSlice();
            voxels.SaveSliceWithTracksAsImage(peakY, "kernel_peak_with_tracks.png", allTracks, photonEnergy);
        }
        private static void SimulatePhoton(WaterTankPhantom voxels, double energy)
        {
            //Launch photons from center of XY at Z=0
            double centerX = (voxels.Nx * voxels.VoxelSizeCm) / 2;
            double centerY = (voxels.Ny * voxels.VoxelSizeCm) / 2;
            Vector3 position = new Vector3((float)centerX, (float)centerY, 0f);
            Vector3 direction = new Vector3(0f, 0f, 0f);

            //Sample distance to interaction
            double mu = PhysicsData.GetTotalLinearAttenuationCoefficient(energy); //in cm-1
            var random = new Random(); var R = random.NextDouble();
            double distance = -(1 / mu) * Math.Log(R);
            Vector3 interactionPoint = position + direction * (float)distance;

            //Check if the photon left the volume
            if (!voxels.BoundaryCheck(interactionPoint.X, interactionPoint.Y, interactionPoint.Z)) 
                return;

            //Select interaction type
            string interaction = PhysicsData.SampleInteractionType(energy);

            //Compute electron parameters and transport it
            var electronFromInteraction = PhysicsData.GetElectronFromInteraction(interaction,energy,direction);
            TransportElectron(voxels,interactionPoint, 
                electronFromInteraction.direction, electronFromInteraction.electronEnergy);
        }
        private static List<List<Vector3>> allTracks = new();
        private static void TransportElectron(
            WaterTankPhantom voxels, 
            Vector3 startPos, Vector3 direction, double energy)
        {
            const double stepSize = 0.05; //in cm
            const double minEnergy = 0.02; //MeV cutoff
            Vector3 position = startPos;
            var currentTrack = new List<Vector3>();
            while (energy > minEnergy && voxels.BoundaryCheck(position.X,position.Y,position.Z))
            {
                //Get stopping power from the CSDA range
                double stoppingPower = CSDArangeHelper.GetStoppingPower(energy);

                //Simulate energy loss for the particle
                double energyLoss = Math.Min(energy, stoppingPower * stepSize);
                energy -= energyLoss;

                //Deposit energy into voxel of interaction
                voxels.DepositEnergy(position.X, position.Y, position.Z, energyLoss);

                //Simulate multiple scattering
                double theta = PhysicsData.SampleMultipleScatteringAngle(energy);
                double phi = PhysicsData.SampleAzimutalAngle();
                Vector3 scatterDirection = PhysicsData.SphericalToCartesian(theta, phi);
                direction = VectorRotationHelper.RotateToDirection(scatterDirection, direction);

                //Take a step
                position += direction * (float)stepSize;
                currentTrack.Add(position);
            }
            allTracks.Add(currentTrack);
        }
    }
}
