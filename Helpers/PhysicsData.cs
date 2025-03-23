using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public static class PhysicsData
    {
        public class CrossSections
        {
            public double Photoelectric;
            public double Compton;
            public double Pair;
            public double Total => Photoelectric + Compton + Pair;
        }
        public static CrossSections GetCrossSections(double energy)
        {
            if (Math.Abs(energy - 2.0) < 0.1)
               return new CrossSections { Photoelectric = 1.063e-6, Compton = 4.901e-2, Pair = 3.908e-4 };
            if (Math.Abs(energy - 6.0) < 0.1)
                return new CrossSections { Photoelectric = 2.483e-7, Compton = 2.455e-2, Pair = 3.155e-3 };
            if (Math.Abs(energy - 10.0) < 0.1)
                return new CrossSections { Photoelectric = 1.386e-7, Compton = 1.710e-2, Pair = 5.090e-3 };

            throw new ArgumentException("Unsupported photon energy");
        }
        public static double GetTotalLinearAttenuationCoefficient(double energy)
        {
            double atomsPerCm3 = 3.34e22; //for water
            var cs = GetCrossSections(energy);
            return cs.Total * 1e-24 * atomsPerCm3;
        }
        public static string SampleInteractionType(double energy)
        {
            var cs = GetCrossSections(energy);
            double totalCs = cs.Total;
            double peProb = cs.Photoelectric/ cs.Total;
            double pairProb = cs.Pair/ cs.Total;
            double comptonProb = cs.Compton/ cs.Total;
            var random = new Random();
            var R = random.NextDouble();
            if (R <= peProb) return "photoelectric";
            if (R > peProb && R <= pairProb) return "pair";
            if (R > pairProb && R <= 1) return "compton";
            else throw new ArgumentException("Check logic for interaction probabilities");
        }
        public static (double electronEnergy, Vector3 direction) GetElectronFromInteraction(
            string interaction, double photonEnergy, Vector3 photonDir)
        {
            double theta = 0;
            double phi = 0;
            Vector3 eDirection = photonDir;
            double eEnergy = 0;

            switch (interaction)
            {
                case "photoelectric":
                    eEnergy = photonEnergy;
                    eDirection = photonDir;
                    break;
                case "compton":
                    phi = SampleComptonAngle(photonEnergy);
                    double tanHalfPhi = Math.Tan(phi / 2.0);
                    double cotTheta = 1 + (photonEnergy / 0.511) * tanHalfPhi;
                    theta = Math.Atan(1 / cotTheta);
                    eEnergy = photonEnergy - ComptonScatteredPhotonEnergy(photonEnergy, phi);
                    eDirection = SphericalToCartesian(theta, SampleAzimutalAngle());
                    break;
                case "pair":
                    double T = photonEnergy - 1.022;
                    eEnergy = SamplePairProductionEnergy(photonEnergy);
                    theta = SamplePairProductionAngle(photonEnergy);
                    phi = SampleAzimutalAngle();
                    eDirection = SphericalToCartesian(theta, phi);
                    break;
            }
            return (eEnergy, eDirection);
        }
        public static double SampleAzimutalAngle()
        {
            var random = new Random();
            var R = random.NextDouble();
            return 2*Math.PI*R;
        }
        public static double ComptonScatteredPhotonEnergy(double energy,double phi)
        {
            double mec2 = 0.511;
            return energy / (1 + (energy / mec2) * (1 - Math.Cos(phi)));
        }
        public static double SamplePairProductionEnergy(double photonEnergy)
        {
            double Tmax = photonEnergy - 1.022;

            const int maxTries = 10000;
            for (int i = 0; i < maxTries; i++)
            {
                var random = new Random();
                var R = random.NextDouble();
                var random1 = new Random();
                var R1 = random1.NextDouble();

                double x = R;
                double fx = Math.Log(1000 * x + 1);
                double y = R1 * Math.Log(1001);
                if (y < fx) return x * Tmax;
            }
            return Tmax / 2;
        }
        public static double SampleComptonAngle(double energy)
        {
            const int maxTries = 10000;
            double mec2 = 0.511;
            for (int i=0;i<maxTries; i++)
            {
                var random = new Random();
                var R = random.NextDouble();
                double phi = R * Math.PI;
                double cosPhi = Math.Cos(phi);
                double ePrime = energy / (1 + (energy / mec2) * (1 - cosPhi));
                double dSigma = Math.Pow(ePrime / energy, 2)*
                    (ePrime/energy + energy/ePrime - Math.Pow(Math.Sin(phi),2));
                if (R < dSigma / 2) return phi;
            }
            return Math.PI / 2;
        }
        public static double SamplePairProductionAngle(double energy)
        {
            double gamma = energy / 0.511 + 1;
            double beta = Math.Sqrt(1 - 1 / (gamma * gamma));
            const int maxTries = 10000;
            for (int i = 0; i < maxTries; i++)
            {
                var random = new Random();
                var R = random.NextDouble();
                double cosTheta = 2 * R - 1;
                double theta = Math.Acos(cosTheta);

                double pdf = 1 / Math.Pow(1 - beta * cosTheta, 2);
                var random1 = new Random();
                var R1 = random1.NextDouble();
                double y = R1 * 4;
                if (y < pdf) return theta;
            }
            return Math.PI / 2;
        }
        public static Vector3 SphericalToCartesian(double theta,double phi)
        {
            double sinTheta = Math.Sin(theta);
            return new Vector3(
                (float)(sinTheta*Math.Cos(phi)),
                (float)(sinTheta * Math.Sin(phi)),
                (float)(Math.Cos(theta)));
        }
        public static double SampleMultipleScatteringAngle (double eMeV)
        {
            double thetaSquared = GetMeanScatteringAngleSquaredForWater(eMeV);

            //Sample from Rayleigh distribution
            var random = new Random(); var R = random.NextDouble();
            return Math.Sqrt(-thetaSquared * Math.Log(R));
        }
        private static double GetMeanScatteringAngleSquaredForWater(double eMeV)
        {
            //Const for water
            double Z = 7.42;
            double A = 14.94;
            double z = 1;
            double m = 0.511;

            //Kinematics
            double Etotal = eMeV + m;
            double gamma = Etotal / m;
            double beta = Math.Sqrt(1.0 - 1.0 / (gamma * gamma));
            double p = Math.Sqrt(eMeV * eMeV + 2 * eMeV * m);

            //Compute chi^2
            double chiSquared = 2.007e-5 * Z * (Z + 1) / (A * p * p * beta * beta);

            //Compute omega and nu
            double Omega = chiSquared * z * z;
            double F = 0.9; //correction factor (empirical)
            double nu = 0.5 * Omega / (1.0 - F);

            //Full formula for thetaSquared
            double logTerm = Math.Log(1 + nu);
            double thetaSquared = (2 * chiSquared) / (1 + F * F) / ((1 + nu) / nu * logTerm - 1);

            return thetaSquared; //in radians
        }

    }
}
