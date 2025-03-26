using MCassignment_PHYS539_PeidusD_2025.Containers;
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
            if (Math.Abs(energy - 2.0) < 0.1)
                return 0.049;
            if (Math.Abs(energy - 6.0) < 0.1)
                return 0.029;
            if (Math.Abs(energy - 10.0) < 0.1)
                return 0.022;

            throw new ArgumentException("Unsupported attenuation coefficient");
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
        public static (double electronEnergy, Vector3 electronDirection, double positronEnergy, Vector3 positronDirection) GetElectronFromInteraction(
            string interaction, double photonEnergy, Vector3 photonDir, 
            StatisticsToExport statistics, double maxKleinNichProb)
        {
            double theta = 0;
            double phi = 0;
            double theta_p = 0;
            double phi_p = 0;
            Vector3 eDirection = photonDir;
            Vector3 pDirection = new Vector3(0,0,0);
            double eEnergy = 0;
            double pEnergy = 0;
            switch (interaction)
            {
                case "photoelectric":
                    eEnergy = photonEnergy;
                    eDirection = photonDir;
                    break;
                case "compton":
                    // Step 1: Sample realistic scattered photon direction
                    Vector3 scatteredPhotonDir = PhysicsData.SampleScatteredComptonPhotonDirection(photonEnergy, photonDir, maxKleinNichProb);

                    // Step 2: Calculate angle between original and scattered direction
                    double anglePhoton = Vector3AngleBetween(photonDir, scatteredPhotonDir);

                    double scatteredEnergy = photonEnergy / (1 + (photonEnergy / 0.511) * (1 - Math.Cos(anglePhoton)));

                    // Step 3: Compute electron energy
                    eEnergy = photonEnergy - scatteredEnergy;

                    // Step 4: Estimate electron angle based on momentum conservation
                    double cotThetaElectron = (1+ photonEnergy / 0.511)*Math.Tan(anglePhoton/2);

                    double thetaElectron = Math.Atan(1/cotThetaElectron);

                    // Step 5: Sample random azimuthal angle for spread
                    double phiElectron = PhysicsData.SampleAzimutalAngle();
                    Vector3 localEDir = PhysicsData.SphericalToCartesian(thetaElectron, phiElectron);
                    eDirection = PhysicsData.RotateToDirection(localEDir, photonDir);
                    statistics.RecordComptonScatteringAngle(anglePhoton);
                    break;
                case "pair":
                    double T = photonEnergy - 1.022;
                    eEnergy = SamplePairProductionEnergy(photonEnergy);
                    pEnergy = T - eEnergy;
                    theta = SamplePairProductionAngle(eEnergy);
                    phi = SampleAzimutalAngle();

                    theta_p = SamplePairProductionAngle(eEnergy);
                    phi_p = SampleAzimutalAngle();

                    eDirection = SphericalToCartesian(theta, phi);
                    pDirection = SphericalToCartesian(theta_p, phi_p);
                    statistics.RecordPairScatteringAngle(theta);
                    break;
            }
            return (eEnergy, eDirection, pEnergy, pDirection);
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
        public static Vector3 SampleScatteredComptonPhotonDirection(double photonEnergyMeV, Vector3 incidentDirection, double maxKleinNichProb)
        {
            const double mec2 = 0.511; // MeV
            Random rng = new Random();
            double theta, phi;
            var electronRadius = Constants.electronRadius_m;

            // Rejection sampling using Klein-Nishina formula
            while (true)
            {
                theta = rng.NextDouble() * Math.PI; // phi_p in [0, pi]

                var hv = photonEnergyMeV;
                var hv_prime = hv / (1+(hv/0.511)*(1-Math.Cos(theta)));

                var thisKlNichCrs = (Math.PI) * (Math.Sin(theta)) * (Math.Pow((hv_prime / hv), 2));
                thisKlNichCrs = thisKlNichCrs * (hv / hv_prime + hv_prime / hv - Math.Pow(Math.Sin(theta), 2));
                thisKlNichCrs = thisKlNichCrs * Math.Pow(electronRadius, 2);

                double maxProb = maxKleinNichProb; // Safe upper bound
                
                if (rng.NextDouble() * maxProb <= thisKlNichCrs)
                    break;
            }

            // Azimuthal angle sampled uniformly
            phi = 2 * Math.PI * rng.NextDouble();

            // Local direction (relative to z-axis)
            Vector3 localDirection = SphericalToCartesian(theta, phi);

            // Rotate local direction to align with incident photon direction
            return RotateToDirection(localDirection, incidentDirection);
        }
        public static Vector3 RotateToDirection(Vector3 local, Vector3 targetDir)
        {
            Vector3 zAxis = new Vector3(0, 0, 1);
            if (Vector3.Distance(Vector3.Normalize(targetDir), zAxis) < 1e-6)
                return local;

            Vector3 axis = Vector3.Cross(zAxis, targetDir);
            float angle = (float)Math.Acos(Vector3.Dot(zAxis, Vector3.Normalize(targetDir)));
            Quaternion rotation = Quaternion.CreateFromAxisAngle(Vector3.Normalize(axis), angle);
            return Vector3.Transform(local, rotation);
        }
        public static double SamplePairProductionAngle(double particleEnergy)
        {
            //Sample scattering angle randomly from 0 to 2PI
            var random = new Random();
            var R = random.NextDouble();
            R = R * 2*Math.PI;
            //double beta = Math.Sqrt(1 - 1 / Math.Pow(((particleEnergy/0.511)+1),2));
            //double scatteringAngleCos = (beta * R - R + 1) / (beta * beta * R - beta * R + 1);
            //double scatteringAngle = Math.Acos(scatteringAngleCos);

            double beta = Math.Sqrt(1 - 1 / Math.Pow((particleEnergy / 0.511) + 1, 2));
            double alarama;
            if (beta > 1 || beta < 0)
                alarama = 0.0;
            double scatAnglCos = (1 / beta) * (1 - (1 / ((1 / (1 - beta)) - (beta * R / (2 * Math.PI)))));
            double scatteringAngle = Math.Acos(scatAnglCos);
            return scatteringAngle;
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
            thetaSquared = 0.01;
            //Sample from Rayleigh distribution via Transformation method
            var random = new Random(); var R = random.NextDouble();
            var multipleScatteringAngle = Math.Sqrt(-thetaSquared * Math.Log(1- R));

            return multipleScatteringAngle;
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
        public static double Vector3AngleBetween(Vector3 a, Vector3 b)
        {
            double dot = Vector3.Dot(Vector3.Normalize(a), Vector3.Normalize(b));
            return Math.Acos(Math.Clamp(dot, -1.0, 1.0));
        }

    }
}
