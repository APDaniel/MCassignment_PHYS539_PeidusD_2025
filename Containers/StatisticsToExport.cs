namespace MCassignment_PHYS539_PeidusD_2025.Containers
{
    public class StatisticsToExport
    {
        public static List<double> scatteringAnglesComptonPhotons = new List<double>();
        public static List<double> scatteringAnglesPairPhotons_e = new List<double>();
        public static List<double> scatteringAnglesPairPhotons_p = new List<double>();
        public static List<double> hardCollisionAngles = new List<double>();
        public static List<double> photonTravelDistance = new List<double>();
        public static double comptonEnergyTransfered;

        //Data for Compton angles
        public List<double> GetPhotoelectricScatteringAngles()
        {
            var test = scatteringAnglesComptonPhotons;
            return test;
        }
        public void RecordComptonScatteringAngle(double angle)
        {
            scatteringAnglesComptonPhotons.Add(angle);
        }

        public void ClearCompton()
        {
            scatteringAnglesComptonPhotons.Clear();
        }
        public double AverageAngleCompton()
        {
            return scatteringAnglesComptonPhotons.Count > 0 ? scatteringAnglesComptonPhotons.Average() : 0.0;
        }

        //Data for electron scattering pair production angles
        public void RecordPairScatteringAngle_e(double angle)
        {
            scatteringAnglesPairPhotons_e.Add(angle);
        }

        public void ClearPair_e()
        {
            scatteringAnglesPairPhotons_e.Clear();
        }
        public double AverageAnglePair_e()
        {
            return scatteringAnglesPairPhotons_e.Count > 0 ? scatteringAnglesPairPhotons_e.Average() : 0.0;
        }
        //Data for electron scattering pair production angles
        public void RecordPairScatteringAngle_p(double angle)
        {
            scatteringAnglesPairPhotons_p.Add(angle);
        }

        public void ClearPair_p()
        {
            scatteringAnglesPairPhotons_p.Clear();
        }
        public double AverageAnglePair_p()
        {
            return scatteringAnglesPairPhotons_p.Count > 0 ? scatteringAnglesPairPhotons_p.Average() : 0.0;
        }

        //Data for electron multiple scattering transport simulation
        public void RecordHardCollisionAngle(double angle)
        {
            hardCollisionAngles.Add(angle);
        }

        public void ClearHardCollisionAngle()
        {
            hardCollisionAngles.Clear();
        }
        public double AverageHardCollisionAngle()
        {
            return hardCollisionAngles.Count > 0 ? hardCollisionAngles.Average() : 0.0;
        }
        

        //Data for electron multiple scattering transport simulation
        public void SavePhotonTravelDistance(double angle)
        {
            photonTravelDistance.Add(angle);
        }

        public void ClearphotonTravel()
        {
            photonTravelDistance.Clear();
        }
        public double AveragePhotonTravel()
        {
            return photonTravelDistance.Count > 0 ? photonTravelDistance.Average() : 0.0;
        }

        //Data for Compton energy transfered to beta-
        public void SaveComptonEnergy(double energy)
        {
            comptonEnergyTransfered = comptonEnergyTransfered + energy; ;
        }

        public void ClearComptonEnergy()
        {
            comptonEnergyTransfered = 0;
        }
    }
}
