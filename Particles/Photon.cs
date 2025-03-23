using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Particles
{
    public class Photon
    {
        public Vector3 Position { get;set; }
        public Vector3 Direction { get;set;}
        public double Energy { get; set; }

        public Photon(Vector3 position, Vector3 direction, double energy)
        {
            Position = position;
            Direction = direction;
            Energy = energy;
        }
    }
}
