using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Helpers
{
    public class VectorRotationHelper
    {
        public static Vector3 RotateToDirection(Vector3 localVec, Vector3 targetDir)
        {
            if (targetDir == Vector3.UnitZ) return localVec;

            Vector3 axis = Vector3.Cross(Vector3.UnitZ, targetDir);
            float angle = (float)Math.Acos(Vector3.Dot(Vector3.UnitZ, Vector3.Normalize(targetDir)));

            if (axis.LengthSquared() < 1e-6) return localVec;

            axis = Vector3.Normalize(axis);
            return RotateAroundAxis(localVec, axis, angle);
        }
        private static Vector3 RotateAroundAxis(Vector3 vec, Vector3 axis, float angle)
        {
            float cos = (float)Math.Cos(angle);
            float sin = (float)Math.Sin(angle);
            return vec * cos + Vector3.Cross(axis, vec) * sin + axis * Vector3.Dot(axis, vec) * (1 - cos);
        }
    }
}
