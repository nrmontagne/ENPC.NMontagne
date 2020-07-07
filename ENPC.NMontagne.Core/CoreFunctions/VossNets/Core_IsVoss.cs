using System;
using System.Collections.Generic;

using ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;
using ENPC.Numerics;
using MathNet.Numerics.Distributions;
using System.Globalization;
using MathNet.Numerics.Providers.LinearAlgebra;

namespace ENPC.NMontagne.Core.CoreFunctions.VossNets
{
    public static class IsVoss
    {
        /// <summary>
        /// Verifies the equality of opposite angles made by edges at vertices with a valences 4. 
        /// </summary>
        /// <param name="mesh"> The mesh to operate on.</param>
        /// <param name="AreFalse"> The list of points which do not respect the equality criteria.</param>
        /// <param name="AreTrue"> The list of points which respect the equality criteria.</param>
        public static void Core_IsGeodesicNet(HeMesh<Point> mesh, out List<Point> AreFalse, out List<Point> AreTrue)
        {
            // Initialization
            AreFalse = new List<Point>();
            AreTrue = new List<Point>();

            foreach (HeVertex<Point> vertex in mesh.GetVertices())
            {
                // There is no criteria on boundary vertices
                if (vertex.IsBoundary()) { continue; }

                List<HeVertex<Point>> neighbours = vertex.NeighbourVertices();
                // The vertex must have four connected edges
                if (neighbours.Count != 4) { throw new ArgumentException("A vertex has less or more than 4 connected edges."); }


                // Defines the vector around the vertex                     
                Vector δF1 = (Vector)(neighbours[0].Position - vertex.Position);
                Vector δF2 = (Vector)(neighbours[1].Position - vertex.Position);
                Vector δF_1 = (Vector)(neighbours[2].Position - vertex.Position);
                Vector δF_2 = (Vector)(neighbours[3].Position - vertex.Position);

                // Computes the angles between the successive vectors.                     
                double α1 = Vector.AngleBetween(δF1, δF2);
                double α2 = Vector.AngleBetween(δF2, δF_1);
                double α3 = Vector.AngleBetween(δF_1, δF_2);
                double α4 = Vector.AngleBetween(δF_2, δF1);

                // Fill results                     
                if (Math.Abs(α1 - α3) < Settings._angularPrecision && Math.Abs(α2 - α4) < Settings._angularPrecision)
                {
                    AreTrue.Add(vertex.Position);
                }
                else { AreFalse.Add(vertex.Position); }
            }
        }

        /// <summary>
        /// Verifies the planarity of the faces. 
        /// </summary>
        /// <param name="mesh"> The mesh to operate on.</param>
        /// <param name="AreFalse"> The list of barycentre of faces which do not respect the equality criteria.</param>
        /// <param name="AreTrue"> The list of barycentre of faces which respect the equality criteria.</param>
        public static void Core_HasPlanarFaces(HeMesh<Point> mesh, out List<Point> AreFalse, out List<Point> AreTrue)
        {
            AreFalse = new List<Point>();
            AreTrue = new List<Point>();

            foreach(HeFace<Point> face in mesh.GetFaces())
            {
                bool isPlanar = true;

                List<HeVertex<Point>> faceVertices = face.FaceVertices();
                int nb_FaceVertex = faceVertices.Count;

                if (nb_FaceVertex != 3)
                {
                    // Create a set of normal vectors from the corss product of edges
                    List<Vector> normals = new List<Vector>();
                    for (int i_Vertex = 0; i_Vertex < nb_FaceVertex - 1; i_Vertex += 2)
                    {
                        int j_Vertex = i_Vertex + 1; int k_Vertex = (i_Vertex + 2) % nb_FaceVertex;
                        Vector dir1 = (Vector)(faceVertices[i_Vertex].Position - faceVertices[j_Vertex].Position);
                        Vector dir2 = (Vector)(faceVertices[k_Vertex].Position - faceVertices[j_Vertex].Position);
                        normals.Add(Vector.CrossProduct(dir1, dir2));
                    }

                    for (int i_Normal = 1; i_Normal < normals.Count; i_Normal++)
                    {
                        if(Vector.CrossProduct(normals[0], normals[i_Normal]).Length() > Settings._absolutePrecision)
                        {
                            isPlanar = false;
                            break;
                        }
                    }
                }

                // Compute the face barycentre
                Point barycenter = new Point();
                for (int i_Vertex = 0; i_Vertex < nb_FaceVertex; i_Vertex++)
                {
                    barycenter += faceVertices[i_Vertex].Position;
                }
                barycenter /= nb_FaceVertex;

                // Store the barycentre
                if (isPlanar) { AreTrue.Add(barycenter); }
                else { AreTrue.Add(barycenter); }

            }
        }
    }
}
