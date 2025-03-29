using System.Drawing;
using System.Drawing.Drawing2D;
using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Containers
{
    public class WaterTankPhantom_ReducedResolution
    {
        private readonly double[,,] doseGrid;
        public int Nx { get; set; }
        public int Ny { get; set; }
        public int Nz { get; set; }
        public double VoxelSizeCm_XY { get; set; }

        public double VoxelSizeCm_Z { get; set; }
        public object Bitmap { get; private set; }

        public WaterTankPhantom_ReducedResolution(int nx, int ny, int nz, double voxelSizeCm_XY, double voxelSizeCm_Z)
        {
            Nx = nx; Ny = ny; Nz = nz; VoxelSizeCm_XY = voxelSizeCm_XY; VoxelSizeCm_Z = voxelSizeCm_Z;
            doseGrid = new double[nx, ny, nz];
        }

        public void DepositEnergy(double x, double y, double z, double energy)
        {
            int ix = (int)(x / VoxelSizeCm_XY);
            int iy = (int)(y / VoxelSizeCm_XY);
            int iz = (int)(z / VoxelSizeCm_Z);

            if (ix >= 0 && ix < Nx &&
                iy >= 0 && iy < Ny &&
                iz >= 0 && iz < Nz)
            {
                doseGrid[ix, iy, iz] += Math.Round(energy, 3);
            }
        }

        public double[,,] GetDoseGrid()
        {
            return doseGrid;
        }
        public void ExportSliceY(int yIndex, string filename, double numParticles, double photonEnergyMeV)
        {
            var normalizedDoseGrid = NormalizeToIncidentEnergy(1, 1);
            using var writer = new StreamWriter(photonEnergyMeV + "_" + filename);

            writer.Write("X\\Z,");
            for (int iz = 0; iz < Nz; iz++)
            {
                double zCoord = iz * VoxelSizeCm_Z;
                writer.Write($"{zCoord:F2};");
            }
            writer.WriteLine();

            for (int ix = 0; ix < Nx; ix++)
            {
                double xCoord = ix * VoxelSizeCm_XY;
                writer.Write($"{xCoord:F2};");

                for (int iz = 0; iz < Nz; iz++)
                {
                    double val = normalizedDoseGrid[ix, yIndex, iz];
                    writer.Write($"{val:E4};");
                }
                writer.WriteLine();
            }
            writer.Close();

            /*string path = "TEST_isolines.csv"; // path to your uploaded CSV

            // === STEP 1: Load CSV ===
            string[] lines = File.ReadAllLines(path);
            ScottPlot.Plot plt = new();
            // Parse Z (column) positions from header
            double[] zCoords = lines[0]
                .Split(';').Skip(1)
                .Select(s => double.Parse(s.Trim()))
                .ToArray();

            int Nz1 = zCoords.Length;
            int Nx1 = lines.Length - 1;

            double[] xCoords = new double[Nx1+1];
            double[,] doseValues = new double[Nx, Nz1];

            for (int ix = 0; ix < Nx1; ix++)
            {
                string[] parts = lines[ix].Split(';');

                xCoords[ix] = double.Parse(parts[0]); // X position
                for (int iz = 0; iz < Nz1; iz++)
                    doseValues[ix, iz] = double.Parse(parts[iz + 1]);
            }

            // === STEP 2: Build Coordinates3d[,] ===
            Coordinates3d[,] cs = new Coordinates3d[Nx1, Nz1];
            for (int ix = 0; ix < Nx1; ix++)
            {
                for (int iz = 0; iz < Nz1; iz++)
                {
                    double x = xCoords[ix];
                    double z = zCoords[iz];
                    double val = doseValues[ix, iz];
                    cs[iz, ix] = new Coordinates3d(x, z, val); // (Z, X, Value)
                }
            }
            double maxVal = cs.Cast<Coordinates3d>().Max(c => c.Z);
            var contour = plt.Add.ContourLines(cs);
            contour.ContourLineLevels = new double[]
            {
                0.2 * maxVal,
                0.4 * maxVal,
                0.6 * maxVal,
                0.8 * maxVal,
                1.0 * maxVal
            };
            var heatmap = plt.Add.Heatmap(cs);
            heatmap.FlipVertically = true;
            heatmap.Colormap = new ScottPlot.Colormaps.MellowRainbow();
            contour.LabelStyle.Bold = true;
            contour.LinePattern = LinePattern.DenselyDashed;
            contour.LineColor = Colors.Black.WithAlpha(.5);
            plt.Axes.TightMargins();
            plt.HideAxesAndGrid();
            plt.SavePng(photonEnergyMeV.ToString() + "MeV" + "_isoLines.png", 1200, 800);*/

            Console.WriteLine($"Exported Y-slice {yIndex} to CSV: {filename}");
        }
        public void ExportFullKernel(string filename, double numPhotons, double photonEnergyMeV)
        {
            using StreamWriter writer = new StreamWriter(filename);
            for (int iy = 0; iy < Ny; iy++)
            {
                double yCoord = iy * VoxelSizeCm_XY;
                var normalizedDoseGrid = NormalizeToIncidentEnergy(numPhotons, photonEnergyMeV);
                writer.WriteLine($"#Slice y={iy} (Y={yCoord:F2} cm)");
                writer.Write("X\\Z");
                for (int iz = 0; iz < Nz; iz++)
                {
                    double zCoord = iz * VoxelSizeCm_Z;
                    writer.Write($";{zCoord:F2},");
                }
                writer.WriteLine();

                for (int ix = 0; ix < Nx; ix++)
                {
                    double xCoord = ix * VoxelSizeCm_XY;
                    writer.Write($"{xCoord:F2}");
                    for (int iz = 0; iz < Nz; iz++)
                    {
                        double val = normalizedDoseGrid[ix, iy, iz];
                        writer.Write($";{val:E4}");
                    }
                    writer.WriteLine();
                }
                writer.WriteLine();
            }
            Console.WriteLine($"Full Kernel exported to: {filename}");
        }
        public double[,,] NormalizeToIncidentEnergy(double numPhotons, double photonEnergyMeV)
        {
            double totalIncidentEnergy = numPhotons * photonEnergyMeV;
            double[,,] normalizedDoseGrid = new double[Nx, Ny, Nz];

            for (int ix = 0; ix < Nx; ix++)
            {
                for (int iy = 0; iy < Ny; iy++)
                {
                    for (int iz = 0; iz < Nz; iz++)
                    {
                        normalizedDoseGrid[ix, iy, iz] = doseGrid[ix, iy, iz] / totalIncidentEnergy;
                    }
                }
            }
            return normalizedDoseGrid;
        }

        public bool BoundaryCheck(double x, double y, double z)
        {
            double sizeX = Nx * VoxelSizeCm_XY;
            double sizeY = Ny * VoxelSizeCm_XY;
            double sizeZ = Nz * VoxelSizeCm_Z;

            return (x >= 0 && x < sizeX &&
                    y >= 0 && y < sizeY &&
                    z >= 0 && z < sizeZ);
        }
        public void SaveSliceWithTracksAsImage(
            int yIndex, string filename,
            List<List<Vector3>> electronTracks,
            double photonEnergy)
        {
            int scale = 1; //test, DNU
            int width = Nx;
            int height = Nz;
            Bitmap bigBmp = new Bitmap(width * scale, height * scale);
            Graphics g = Graphics.FromImage(bigBmp);
            g.SmoothingMode = SmoothingMode.HighQuality;
            g.Clear(System.Drawing.Color.LightBlue);
            //Draw tracks
            Pen trackPen = new Pen(System.Drawing.Color.Black, 1);

            foreach (var track in electronTracks)
                for (int i = 1; i < track.Count; i++)
                {
                    var p1 = track[i - 1];
                    var p2 = track[i];

                    //Plot only if both points are in the same Y-slice
                    int iy1 = (int)(p1.Y / VoxelSizeCm_XY);
                    int iy2 = (int)(p2.Y / VoxelSizeCm_XY);
                    if (iy1 == yIndex && iy2 == yIndex)
                    {
                        int x1 = (int)(p1.X / VoxelSizeCm_XY);
                        int z1 = (int)(p1.Z / VoxelSizeCm_Z);
                        int x2 = (int)(p2.X / VoxelSizeCm_XY);
                        int z2 = (int)(p2.Z / VoxelSizeCm_Z);

                        g.DrawLine(trackPen, x1, z1, x2, z2);
                    }
                }
            Bitmap finalBitmap = new Bitmap(Nx, Nz);
            Graphics gFinal = Graphics.FromImage(finalBitmap);
            gFinal.DrawImage(bigBmp, 0, 0, Nx, Nz);
            finalBitmap.Save($"{photonEnergy}MeV_" + filename, System.Drawing.Imaging.ImageFormat.Png);
            Console.WriteLine($"Saved slice+electron tracks: {filename}");
        }

        public int FindPeakYSlice()
        {
            double max = 0;
            int bestY = Ny / 2;
            for (int iy = 0; iy < Ny; iy++)
            {
                double sum = 0;
                for (int ix = 0; ix < Nx; ix++)
                    for (int iz = 0; iz < Nz; iz++)
                        sum += doseGrid[ix, iy, iz];
                if (sum > max)
                {
                    max = sum;
                    bestY = iy;
                }
            }
            return bestY;
        }
        public double[,] GetSlice2D_Y(int yIndex)
        {
            double[,] slice = new double[Nx, Nz];

            for (int ix = 0; ix < Nx; ix++)
                for (int iz = 0; iz < Nz; iz++)
                    slice[ix, iz] = doseGrid[ix, yIndex, iz];
            return slice;
        }
    }
}
