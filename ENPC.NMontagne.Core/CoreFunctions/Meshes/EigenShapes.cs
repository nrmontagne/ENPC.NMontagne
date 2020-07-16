using System;
using System.Linq;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;
using ENPC.Numerics;

using Acc = Accord.Math;


namespace ENPC.NMontagne.Core.CoreFunctions.Meshes
{
    /// <summary>
    /// Class containing methods to compute the eigenshapes from the gauss map of a mesh.
    /// </summary>
    public static class Eigenshapes
    {
        /// <summary>
        /// Computes the eigenshapes of a net from its Gauss map.
        /// </summary>
        /// <param name="gaussMap"> The Gauss map of the net.</param>
        /// <param name="net"> The net whose eigenshapes are being computed. </param>
        /// <param name="modes"> The eigenshapes of the net.</param>
        public static void Core_FromGaussMap(HeMesh<Euc.Point> gaussMap, out HeMesh<Euc.Point> net, out Dictionary<int, Euc.Vector[]> modes)
        {
            // On the output net
            net = new HeMesh<Euc.Point>();

            int nbVertex;       // Equal to the number of faces in the Gauss image.
            int nbEdge;         // Equal to the number of internal edges in the Gauss image.
            int nbFace;         // Equal to the number of internal vertices in the Gauss image.

            List<int> e_VtoG = new List<int>();     // Convert the net edge indices into Gauss Halfedge indices
            List<int> f_VtoG = new List<int>();     // Convert the net face indices into Gauss vertex indices

            #region Initialisation

            // Number of vertices in the net
            nbVertex = gaussMap.FaceCount;

            // Number of edges in the net
            for (int i = 0; i < gaussMap.EdgeCount / 2; i++)
            {
                if (gaussMap.GetEdge(2 * i).IsBoundary() || gaussMap.GetEdge(2 * i).PairEdge.IsBoundary()) { continue; }
                e_VtoG.Add(2 * i);
            }
            nbEdge = e_VtoG.Count;

            //Number of faces in the net
            for (int i = 0; i < gaussMap.VertexCount; i++)
            {
                if (gaussMap.GetVertex(i).IsBoundary()) { continue; }
                f_VtoG.Add(i);
            }
            nbFace = f_VtoG.Count;

            #endregion

            #region Edge Directions 

            Euc.Vector[] dir_Edges = new Euc.Vector[nbEdge];
            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = Euc.Vector.CrossProduct((Euc.Vector)gaussMap.GetEdge(e_VtoG[i_Edge]).EndVertex.Position,
                    (Euc.Vector)gaussMap.GetEdge(e_VtoG[i_Edge]).StartVertex.Position);

                dir_Edges[i_Edge] = dir;
            }

            #endregion

            #region Edges from Vertices : EfromV

            // With this matrix we have E = EfromV * V.

            double[,] EfromV = Acc.Matrix.Create(3 * nbEdge, 3 * nbVertex, 0.0);

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                int i_StartVertex = gaussMap.GetEdge(e_VtoG[i_Edge]).AdjacentFace.Index;
                int i_EndVertex = gaussMap.GetEdge(e_VtoG[i_Edge]).PairEdge.AdjacentFace.Index;

                // "X component of Edge" = "X component of EndVertex" - "X component of StartVertex" 
                EfromV[(3 * i_Edge), (3 * i_EndVertex)] = 1.0; EfromV[(3 * i_Edge), (3 * i_StartVertex)] = -1.0;
                // "Y component of Edge" = "Y component of EndVertex" - "Y component of StartVertex" 
                EfromV[(3 * i_Edge) + 1, (3 * i_EndVertex) + 1] = 1.0; EfromV[(3 * i_Edge) + 1, (3 * i_StartVertex) + 1] = -1.0;
                // "Z component of Edge" = "Z component of EndVertex" - "Z component of StartVertex" 
                EfromV[(3 * i_Edge) + 2, (3 * i_EndVertex) + 2] = 1.0; EfromV[(3 * i_Edge) + 2, (3 * i_StartVertex) + 2] = -1.0;
            }

            #endregion

            #region Cross Product with Direction : C

            // With this matrix we want : C * EfromV * V = 0. It expresses to goals :
            // 1. The edge "EfromV * V" of the net must be parallel to the direction "dir" computed thanks to the Gauss map.
            // 2. The length of the edges "EfromV * V" around a face must be compatible, ie form a closed face. 
            double[,] C = Acc.Matrix.Create(3 * nbEdge, 3 * nbEdge, 0.0);

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                Euc.Vector dir = dir_Edges[i_Edge];
                // For the X component of the cross product of dir and edge E[i]
                C[(3 * i_Edge), (3 * i_Edge)] = 0.0; C[(3 * i_Edge), (3 * i_Edge) + 1] = -dir.Z; C[(3 * i_Edge), (3 * i_Edge) + 2] = dir.Y;
                // For the Y component of the cross product of dir and edge E[i]
                C[(3 * i_Edge) + 1, (3 * i_Edge)] = dir.Z; C[(3 * i_Edge) + 1, (3 * i_Edge) + 1] = 0.0; C[(3 * i_Edge) + 1, (3 * i_Edge) + 2] = -dir.X;
                // For the Z component of the cross product of dir and edge E[i]
                C[(3 * i_Edge) + 2, (3 * i_Edge)] = -dir.Y; C[(3 * i_Edge) + 2, (3 * i_Edge) + 1] = dir.X; C[(3 * i_Edge) + 2, (3 * i_Edge) + 2] = 0.0;
            }

            #endregion

            #region Null space : N

            // To filter global translations in the modes : We choose to have the vertex "0" at (0,0,0)
            // This results in removing the choice of the first vertex, thus we adapt the matrix EfromV

            double[,] EfromV_f = new double[EfromV.GetLength(0), EfromV.GetLength(1) - 3];
            for (int row = 0; row < EfromV.GetLength(0); row++)
            {
                for (int column = 0; column < EfromV.GetLength(1) - 3; column++)
                {
                    EfromV_f[row, column] = EfromV[row, column + 3];
                }
            }

            // Compute the linear space where V (the oordinate column-vector) can be found.
            // It corresponds to the kernel K of C * EfromV_f).
            double[,] B = Acc.Matrix.Dot(C, EfromV_f);

            // Get the null space of the matrix "prod"
            var svd = new Accord.Math.Decompositions.SingularValueDecomposition(B);

            var svd_V = svd.RightSingularVectors;

            var svd_Threshold = Acc.Matrix.GetLength(B).Max() * Acc.Constants.DoubleEpsilon;
            var E = Acc.Matrix.Find(svd.Diagonal, x => (Double)Math.Abs(x) < svd_Threshold);


            double[,] Kern_f = Acc.Matrix.GetColumns(svd_V, E);



            // Add again the vertices
            double[,] N = Acc.Matrix.Create(Kern_f.GetLength(0) + 3, Kern_f.GetLength(1), 0.0);
            for (int row = 0; row < Kern_f.GetLength(0); row++)
            {
                for (int column = 0; column < Kern_f.GetLength(1); column++)
                {
                    N[row + 3, column] = Kern_f[row, column];
                }
            }

            #endregion

            #region Laplacian : L

            double[,] L = Acc.Matrix.Create(3 * nbVertex, 3 * nbVertex, 0.0);

            for (int i_Vertex = 0; i_Vertex < nbVertex; i_Vertex++)
            {
                List<int> i_Neighbours = new List<int>();
                foreach (HeFace<Euc.Point> face in gaussMap.GetFace(i_Vertex).AdjacentFaces())
                {
                    if (face is null) { continue; }
                    else
                    {
                        i_Neighbours.Add(face.Index);
                    }
                }

                // For the current vertex
                L[(3 * i_Vertex), (3 * i_Vertex)] = 1.0; L[(3 * i_Vertex) + 1, (3 * i_Vertex) + 1] = 1.0; L[(3 * i_Vertex) + 2, (3 * i_Vertex) + 2] = 1.0;

                // For each neighbour vertices
                double weigth = -1.0 / ((double)i_Neighbours.Count);

                foreach (int i_Neighbour in i_Neighbours)
                {
                    L[(3 * i_Vertex), (3 * i_Neighbour)] = weigth; L[(3 * i_Vertex) + 1, (3 * i_Neighbour) + 1] = weigth; L[(3 * i_Vertex) + 2, (3 * i_Neighbour) + 2] = weigth;
                }
            }

            #endregion

            #region Matrix : P

            double[,] Nt = Acc.Matrix.Transpose(N);
            double[,] LN = Acc.Matrix.Dot(L, N);
            double[,] NtLN = Acc.Matrix.Dot(Nt, LN);

            Acc.Decompositions.EigenvalueDecomposition eigenDec = new Acc.Decompositions.EigenvalueDecomposition(NtLN, false, true);

            #endregion

            #region Generate HeMesh

            net = new HeMesh<Euc.Point>();
            for (int i_Vertex = 0; i_Vertex < nbVertex; i_Vertex++)
            {
                net.AddVertex(new Euc.Point(0.0, 0.0, 0.0));
            }
            for (int i_Face = 0; i_Face < nbFace; i_Face++)
            {
                List<HeVertex<Euc.Point>> faceVertices = new List<HeVertex<Euc.Point>>();
                foreach (HeEdge<Euc.Point> halfedge in gaussMap.GetVertex(f_VtoG[i_Face]).OutgoingEdges())
                {
                    faceVertices.Add(net.GetVertex(halfedge.AdjacentFace.Index));
                }
                net.AddFace(faceVertices);
            }

            #endregion

            #region Saving Modes
            double[,] eigenVectors = Acc.Matrix.Dot(N, eigenDec.Eigenvectors);

            modes = new Dictionary<int, Euc.Vector[]>();
            for (int i_modes = 0; i_modes < eigenVectors.GetLength(1); i_modes++)
            {
                modes.Add(i_modes, new Euc.Vector[nbVertex]);
            }

            for (int i_modes = 0; i_modes < eigenVectors.GetLength(1); i_modes++)
            {
                for (int i_Vertex = 1; i_Vertex < nbVertex; i_Vertex++)
                {
                    modes[i_modes][i_Vertex] = new Euc.Vector(eigenVectors[(3 * i_Vertex), i_modes], eigenVectors[(3 * i_Vertex) + 1, i_modes], eigenVectors[(3 * i_Vertex) + 2, i_modes]);
                }
            }

            #endregion
        }

        /// <summary>
        /// Computes the eigenshapes of a net.
        /// </summary>
        /// <param name="net"> The net whose eigenshapes are being computed. </param>
        /// <param name="modes"> The eigenshapes of the net.</param>
        public static void Core_FromNet(HeMesh<Euc.Point> net, out Dictionary<int, Euc.Vector[]> modes)
        {
            int nbVertex = net.VertexCount;       // Equal to the number of faces in the Gauss image.
            int nbEdge = net.EdgeCount;         // Equal to the number of internal edges in the Gauss image.

            #region Edges from Vertices : EfromV

            // With this matrix we have E = EfromV * V.

            double[,] EfromV = Acc.Matrix.Create(3 * nbEdge, 3 * nbVertex, 0.0);

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                HeEdge<Euc.Point> edge = net.GetEdge(i_Edge);
                int i_StartVertex = edge.StartVertex.Index;
                int i_EndVertex = edge.EndVertex.Index;

                // "X component of Edge" = "X component of EndVertex" - "X component of StartVertex" 
                EfromV[(3 * i_Edge), (3 * i_EndVertex)] = 1.0; EfromV[(3 * i_Edge), (3 * i_StartVertex)] = -1.0;
                // "Y component of Edge" = "Y component of EndVertex" - "Y component of StartVertex" 
                EfromV[(3 * i_Edge) + 1, (3 * i_EndVertex) + 1] = 1.0; EfromV[(3 * i_Edge) + 1, (3 * i_StartVertex) + 1] = -1.0;
                // "Z component of Edge" = "Z component of EndVertex" - "Z component of StartVertex" 
                EfromV[(3 * i_Edge) + 2, (3 * i_EndVertex) + 2] = 1.0; EfromV[(3 * i_Edge) + 2, (3 * i_StartVertex) + 2] = -1.0;
            }

            #endregion

            #region Cross Product with Direction : C

            // With this matrix we want : C * EfromV * V = 0. It expresses to goals :
            // 1. The edge "EfromV * V" of the net must be parallel to the direction "dir" computed thanks to the Gauss map.
            // 2. The length of the edges "EfromV * V" around a face must be compatible, ie form a closed face. 
            double[,] C = Acc.Matrix.Create(3 * nbEdge, 3 * nbEdge, 0.0);

            for (int i_Edge = 0; i_Edge < nbEdge; i_Edge++)
            {
                HeEdge<Euc.Point> edge = net.GetEdge(i_Edge);
                Euc.Vector dir = (Euc.Vector) (edge.EndVertex.Position = edge.StartVertex.Position);
                // For the X component of the cross product of dir and edge E[i]
                C[(3 * i_Edge), (3 * i_Edge)] = 0.0; C[(3 * i_Edge), (3 * i_Edge) + 1] = -dir.Z; C[(3 * i_Edge), (3 * i_Edge) + 2] = dir.Y;
                // For the Y component of the cross product of dir and edge E[i]
                C[(3 * i_Edge) + 1, (3 * i_Edge)] = dir.Z; C[(3 * i_Edge) + 1, (3 * i_Edge) + 1] = 0.0; C[(3 * i_Edge) + 1, (3 * i_Edge) + 2] = -dir.X;
                // For the Z component of the cross product of dir and edge E[i]
                C[(3 * i_Edge) + 2, (3 * i_Edge)] = -dir.Y; C[(3 * i_Edge) + 2, (3 * i_Edge) + 1] = dir.X; C[(3 * i_Edge) + 2, (3 * i_Edge) + 2] = 0.0;
            }

            #endregion

            #region Null space : N

            // To filter global translations in the modes : We choose to have the vertex "0" at (0,0,0)
            // This results in removing the choice of the first vertex, thus we adapt the matrix EfromV

            double[,] EfromV_f = new double[EfromV.GetLength(0), EfromV.GetLength(1) - 3];
            for (int row = 0; row < EfromV.GetLength(0); row++)
            {
                for (int column = 0; column < EfromV.GetLength(1) - 3; column++)
                {
                    EfromV_f[row, column] = EfromV[row, column + 3];
                }
            }

            // Compute the linear space where V (the oordinate column-vector) can be found.
            // It corresponds to the kernel K of C * EfromV_f).
            double[,] B = Acc.Matrix.Dot(C, EfromV_f);

            // Get the null space of the matrix "prod"
            var svd = new Accord.Math.Decompositions.SingularValueDecomposition(B);

            var svd_V = svd.RightSingularVectors;

            var svd_Threshold = Acc.Matrix.GetLength(B).Max() * Acc.Constants.DoubleEpsilon;
            var E = Acc.Matrix.Find(svd.Diagonal, x => (Double)Math.Abs(x) < svd_Threshold);


            double[,] Kern_f = Acc.Matrix.GetColumns(svd_V, E);



            // Add again the vertices
            double[,] N = Acc.Matrix.Create(Kern_f.GetLength(0) + 3, Kern_f.GetLength(1), 0.0);
            for (int row = 0; row < Kern_f.GetLength(0); row++)
            {
                for (int column = 0; column < Kern_f.GetLength(1); column++)
                {
                    N[row + 3, column] = Kern_f[row, column];
                }
            }

            #endregion

            #region Laplacian : L

            double[,] L = Acc.Matrix.Create(3 * nbVertex, 3 * nbVertex, 0.0);

            for (int i_Vertex = 0; i_Vertex < nbVertex; i_Vertex++)
            {
                List<int> i_Neighbours = new List<int>();
                foreach (HeVertex<Euc.Point> vertex in net.GetVertex(i_Vertex).NeighbourVertices())
                {
                    i_Neighbours.Add(vertex.Index);
                }

                // For the current vertex
                L[(3 * i_Vertex), (3 * i_Vertex)] = 1.0; L[(3 * i_Vertex) + 1, (3 * i_Vertex) + 1] = 1.0; L[(3 * i_Vertex) + 2, (3 * i_Vertex) + 2] = 1.0;

                // For each neighbour vertices
                double weigth = -1.0 / ((double)i_Neighbours.Count);

                foreach (int i_Neighbour in i_Neighbours)
                {
                    L[(3 * i_Vertex), (3 * i_Neighbour)] = weigth; L[(3 * i_Vertex) + 1, (3 * i_Neighbour) + 1] = weigth; L[(3 * i_Vertex) + 2, (3 * i_Neighbour) + 2] = weigth;
                }
            }

            #endregion

            #region Matrix : P

            double[,] Nt = Acc.Matrix.Transpose(N);
            double[,] LN = Acc.Matrix.Dot(L, N);
            double[,] NtLN = Acc.Matrix.Dot(Nt, LN);

            Acc.Decompositions.EigenvalueDecomposition eigenDec = new Acc.Decompositions.EigenvalueDecomposition(NtLN, false, true);

            #endregion

            #region Saving Modes
            double[,] eigenVectors = Acc.Matrix.Dot(N, eigenDec.Eigenvectors);

            modes = new Dictionary<int, Euc.Vector[]>();
            for (int i_modes = 0; i_modes < eigenVectors.GetLength(1); i_modes++)
            {
                modes.Add(i_modes, new Euc.Vector[nbVertex]);
            }

            for (int i_modes = 0; i_modes < eigenVectors.GetLength(1); i_modes++)
            {
                for (int i_Vertex = 1; i_Vertex < nbVertex; i_Vertex++)
                {
                    modes[i_modes][i_Vertex] = new Euc.Vector(eigenVectors[(3 * i_Vertex), i_modes], eigenVectors[(3 * i_Vertex) + 1, i_modes], eigenVectors[(3 * i_Vertex) + 2, i_modes]);
                }
            }

            #endregion
        }
    }
}
