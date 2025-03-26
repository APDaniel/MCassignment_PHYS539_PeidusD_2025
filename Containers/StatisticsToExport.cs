namespace MCassignment_PHYS539_PeidusD_2025.Containers
{
    public class StatisticsToExport
    {
        private List<double> scatteringAnglesComptonPhotons = new List<double>();
        private List<double> scatteringAnglesPairPhotons = new List<double>();

        //Data for Compton angles
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
        public void RecordPairScatteringAngle(double angle)
        {
            scatteringAnglesPairPhotons.Add(angle);
        }

        public void ClearPair()
        {
            scatteringAnglesPairPhotons.Clear();
        }
        public double AverageAnglePair()
        {
            return scatteringAnglesPairPhotons.Count > 0 ? scatteringAnglesPairPhotons.Average() : 0.0;
        }
    }
}
