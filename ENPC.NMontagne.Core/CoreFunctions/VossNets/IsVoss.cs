using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

namespace ENPC.NMontagne.Core.CoreFunctions.VossNets
{
    /// <summary>
    /// Class containing methods to verify that the net is a Voss net.
    /// </summary>
    public static class IsVoss
    {
        /// <summary>
        /// Verifies the equality of opposite angles made by edges at vertices with a valences 4. 
        /// </summary>
        /// <param name="mesh"> The mesh to operate on.</param>
        /// <param name="AreFalse"> The list of points which do not respect the equality criteria.</param>
        /// <param name="AreTrue"> The list of points which respect the equality criteria.</param>
        public static void Core_IsGeodesicNet(HeMesh<Euc.Point> mesh, out List<Euc.Point> AreFalse, out List<Euc.Point> AreTrue)
        {
            // Initialization
            AreFalse = new List<Euc.Point>();
            AreTrue = new List<Euc.Point>();

            foreach (HeVertex<Euc.Point> vertex in mesh.GetVertices())
            {
                // There is no criteria on boundary vertices
                if (vertex.IsBoundary()) { continue; }

                List<HeVertex<Euc.Point>> neighbours = vertex.NeighbourVertices();
                // The vertex must have four connected edges
                if (neighbours.Count != 4) { throw new ArgumentException("A vertex has less or more than 4 connected edges."); }


                // Defines the vector around the vertex                     
                Euc.Vector δF1 = (Euc.Vector)(neighbours[0].Position - vertex.Position);
                Euc.Vector δF2 = (Euc.Vector)(neighbours[1].Position - vertex.Position);
                Euc.Vector δF_1 = (Euc.Vector)(neighbours[2].Position - vertex.Position);
                Euc.Vector δF_2 = (Euc.Vector)(neighbours[3].Position - vertex.Position);

                // Computes the angles between the successive vectors.                     
                double α1 = Euc.Vector.AngleBetween(δF1, δF2);
                double α2 = Euc.Vector.AngleBetween(δF2, δF_1);
                double α3 = Euc.Vector.AngleBetween(δF_1, δF_2);
                double α4 = Euc.Vector.AngleBetween(δF_2, δF1);

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
        public static void Core_HasPlanarFaces(HeMesh<Euc.Point> mesh, out List<Euc.Point> AreFalse, out List<Euc.Point> AreTrue)
        {
            AreFalse = new List<Euc.Point>();
            AreTrue = new List<Euc.Point>();

            foreach(HeFace<Euc.Point> face in mesh.GetFaces())
            {
                bool isPlanar = true;

                List<HeVertex<Euc.Point>> faceVertices = face.FaceVertices();
                int nb_FaceVertex = faceVertices.Count;

                if (nb_FaceVertex != 3)
                {
                    // Create a set of normal vectors from the corss product of edges
                    List<Euc.Vector> normals = new List<Euc.Vector>();
                    for (int i_Vertex = 0; i_Vertex < nb_FaceVertex - 1; i_Vertex += 2)
                    {
                        int j_Vertex = i_Vertex + 1; int k_Vertex = (i_Vertex + 2) % nb_FaceVertex;
                        Euc.Vector dir1 = (Euc.Vector)(faceVertices[i_Vertex].Position - faceVertices[j_Vertex].Position);
                        Euc.Vector dir2 = (Euc.Vector)(faceVertices[k_Vertex].Position - faceVertices[j_Vertex].Position);
                        normals.Add(Euc.Vector.CrossProduct(dir1, dir2));
                    }

                    for (int i_Normal = 1; i_Normal < normals.Count; i_Normal++)
                    {
                        if(Euc.Vector.CrossProduct(normals[0], normals[i_Normal]).Length() > Settings._absolutePrecision)
                        {
                            isPlanar = false;
                            break;
                        }
                    }
                }

                // Compute the face barycentre
                Euc.Point barycenter = new Euc.Point();
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
