using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Containers
{
    public class WaterTankPhantom
    {
        private readonly double[,,] doseGrid;
        public int Nx { get; set; }
        public int Ny { get; set; }
        public int Nz { get; set; }
        public double VoxelSizeCm { get; set; }
        public object Bitmap { get; private set; }

        public WaterTankPhantom(int nx, int ny, int nz, double voxelSizeCm)
        {
            Nx = nx; Ny = ny; Nz = nz; VoxelSizeCm = voxelSizeCm;
            doseGrid = new double[nx, ny, nz];
        }

        public void DepositEnergy(double x, double y, double z, double energy)
        {
            int ix = (int)(x / VoxelSizeCm);
            int iy = (int)(y / VoxelSizeCm);
            int iz = (int)(z / VoxelSizeCm);

            if (ix>=0 && ix<Nx &&
                iy>=0 && iy<Ny &&
                iz>=0 && iz < Nz)
            {
                doseGrid[ix, iy, iz] += energy;
            }
        }

        public double[,,] GetDoseGrid()
        {
            return doseGrid;
        }
        public void ExportSliceY(int yIndex, string filename, double numPhotons, double photonEnergyMeV)
        {
            var normalizedDoseGrid = NormalizeToIncidentEnergy(numPhotons, photonEnergyMeV);
            using var writer = new StreamWriter(filename);

            writer.Write("X\\Z,");
            for (int iz = 0; iz < Nz; iz++)
            {
                double zCoord = iz * VoxelSizeCm;
                writer.Write($"{zCoord:F2},");
            }
            writer.WriteLine();

            for (int ix = 0; ix < Nx; ix++)
            {
                double xCoord = ix * VoxelSizeCm;
                writer.Write($"{xCoord:F2},");

                for (int iz = 0; iz < Nz; iz++)
                {
                    double val = normalizedDoseGrid[ix, yIndex, iz];
                    writer.Write($"{val:E4},");
                }
                writer.WriteLine();
            }
            Console.WriteLine($"Exported Y-slice {yIndex} to CSV: {filename}");
        }
        public void ExportFullKernel(string filename, double numPhotons, double photonEnergyMeV)
        {
            using StreamWriter writer = new StreamWriter(filename);
            for (int iy = 0; iy < Ny; iy++)
            {
                double yCoord = iy * VoxelSizeCm;
                var normalizedDoseGrid = NormalizeToIncidentEnergy(numPhotons, photonEnergyMeV);
                writer.WriteLine($"#Slice y={iy} (Y={yCoord:F2} cm)");
                writer.Write("X\\Z");
                for (int iz = 0; iz < Nz; iz++)
                {
                    double zCoord = iz * VoxelSizeCm;
                    writer.Write($";{zCoord:F2},");
                }
                writer.WriteLine();

                for (int ix = 0; ix < Nx; ix++)
                {
                    double xCoord = ix * VoxelSizeCm;
                    writer.Write($"{xCoord:F2}");
                    for (int iz=0; iz < Nz; iz++)
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
        public double [,,] NormalizeToIncidentEnergy(double numPhotons, double photonEnergyMeV)
        {
            double totalIncidentEnergy = numPhotons * photonEnergyMeV;
            double[,,] normalizedDoseGrid=new double[Nx,Ny,Nz];

            for (int ix=0;ix<Nx; ix++)
            {
                for (int iy = 0; iy < Ny; iy++)
                {
                    for (int iz = 0; iz < Nz; iz++)
                    {
                        normalizedDoseGrid[ix,iy,iz]=doseGrid[ix, iy, iz] /totalIncidentEnergy;
                    }
                }
            }
            return normalizedDoseGrid;
        }

        public bool BoundaryCheck(double x,double y, double z)
        {
            double sizeX = Nx * VoxelSizeCm;
            double sizeY = Ny * VoxelSizeCm;
            double sizeZ = Nz * VoxelSizeCm;

            return (x >= 0 && x < sizeX &&
                    y >= 0 && y < sizeY &&
                    z >= 0 && z < sizeZ);
        }
        public void SaveSliceWithTracksAsImage(
            int yIndex,string filename,
            List<List<Vector3>> electronTracks,
            double photonEnergy)
        {
            int scale = 1;
            int width = Nx;
            int height = Nz;
            Bitmap bigBmp = new Bitmap(width*scale, height*scale);
            Graphics g = Graphics.FromImage(bigBmp);
            g.SmoothingMode = SmoothingMode.HighQuality;
            g.Clear(Color.White);
            //Draw tracks
            Pen trackPen = new Pen(Color.Black,1);

            foreach (var track in electronTracks)
                for (int i = 1; i < track.Count; i++)
                {
                    var p1 = track[i - 1];
                    var p2 = track[i];

                    //Plot only if both points are in the same Y-slice
                    int iy1 = (int)(p1.Y / VoxelSizeCm);
                    int iy2 = (int)(p2.Y / VoxelSizeCm);
                    if (iy1 == yIndex && iy2 == yIndex)
                    {
                        int x1 = (int)(p1.X / VoxelSizeCm);
                        int z1 = (int)(p1.Z / VoxelSizeCm);
                        int x2 = (int)(p2.X / VoxelSizeCm);
                        int z2 = (int)(p2.Z / VoxelSizeCm);

                        g.DrawLine(trackPen, x1, z1, x2, z2);
                    }
                }
            Bitmap finalBitmap = new Bitmap(Nx, Nz);
            Graphics gFinal = Graphics.FromImage(finalBitmap);
            gFinal.DrawImage(bigBmp,0,0,Nx,Nz);
            finalBitmap.Save($"{photonEnergy}MeV_"+filename, ImageFormat.Png);
            Console.WriteLine($"Saved slice+electron tracks: {filename}");
        }
        public void SaveKernelWithIsolines(int iIndex,string filename)
        {

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
    }
}
