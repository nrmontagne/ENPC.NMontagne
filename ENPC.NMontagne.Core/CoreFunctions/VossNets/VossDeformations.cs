using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.Geometry.Euclidean.Extensions.PolyhedralMesh;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;


namespace ENPC.NMontagne.Core.CoreFunctions.VossNets
{
    /// <summary>
    /// Class containing methods to generate Voss nets from its Gauss map using linear optimisation.
    /// </summary>
    public static class VossDeformations
    {
        /// <summary>
        /// Computes the isometric deformation of a Voss net.
        /// </summary>
        /// <param name="initialMesh"> The Voss net to deform.</param>
        /// <param name="factor"> The factor for the deformation.</param>
        /// <param name="deformedMesh"> The deformed Voss net.</param>
        public static void Core_Isometric(HeMesh<Euc.Point> initialMesh, double factor, out HeMesh<Euc.Point> deformedMesh)
        {
            #region Instanciation

            // Prepare meshes for the isometry.
            initialMesh.CleanMesh();
            deformedMesh = (HeMesh<Euc.Point>)initialMesh.Clone();

            // Allows to evaluate whether the vertex or the face as already been modified or not.
            // In this way they will not be modified twice.
            bool[] isVertexModified = new bool[initialMesh.VertexCount];
            bool[] isFaceModified = new bool[initialMesh.FaceCount];

            // Dictionnary fir the iterative modification of the mesh
            // Keys: Index of the face to modify; Value: index of the previous face (which led to the key face);
            Dictionary<int, int> i_CurrFace_Warp = new Dictionary<int, int>();
            Dictionary<int, int> i_NextFace_Warp = new Dictionary<int, int>();
            Dictionary<int, int> i_CurrFace_Weft = new Dictionary<int, int>();
            Dictionary<int, int> i_NextFace_Weft = new Dictionary<int, int>();

            #endregion

            #region Initialize Iteration

            {
                HeFace<Euc.Point> face = deformedMesh.GetFace(0);
                List<HeVertex<Euc.Point>> faceVertices = face.FaceVertices();

                // Set vertices as modified 
                foreach (HeVertex<Euc.Point> faceVertex in faceVertices)
                {
                    isVertexModified[faceVertex.Index] = true;
                }

                // Set face as built
                isFaceModified[face.Index] = true;

                // Set next iteration


                List<HeEdge<Euc.Point>> edges = face.FaceEdges();
                for(int i_Edge = 0; i_Edge < edges.Count; i_Edge++)
                {
                    HeEdge<Euc.Point> edge = edges[i_Edge];
                    if (edge.PairEdge.AdjacentFace == null) { continue; }

                    // Get the index of the face to build next
                    int i_NextFace = edge.PairEdge.AdjacentFace.Index;

                    // Set the next face as modified (eventhough it is not yet appended)
                    isFaceModified[i_NextFace] = true;

                    // Assign face to the Warp or Weft group.
                    if (i_Edge % 2 == 0) { i_NextFace_Warp.Add(i_NextFace, face.Index); }
                    else { i_NextFace_Weft.Add(i_NextFace, face.Index); }
                }
            }

            #endregion

            #region Iteration

            int securityCount = 0;
            while( (i_NextFace_Warp.Count != 0 || i_NextFace_Weft.Count != 0) & securityCount < 2 * initialMesh.FaceCount)
            {
                #region Iteration: Prepare 

                securityCount++;

                // Offset informations
                i_CurrFace_Warp = i_NextFace_Warp;
                i_NextFace_Warp = new Dictionary<int, int>();

                i_CurrFace_Weft = i_NextFace_Weft;
                i_NextFace_Weft = new Dictionary<int, int>();

                #endregion

                #region Iteration : Modify Vertices and Faces

                foreach (int i_Face in i_CurrFace_Warp.Keys)
                {
                    // Move the face vertices.
                    List<HeEdge<Euc.Point>> edges = MoveVertices(deformedMesh, i_CurrFace_Warp[i_Face], i_Face, factor);

                    int nb_Edges = edges.Count;
                    for (int i_Edge = 1; i_Edge < nb_Edges; i_Edge++)
                    {
                        HeFace<Euc.Point> adjacentFace = edges[i_Edge].PairEdge.AdjacentFace;
                        if (adjacentFace == null || isFaceModified[adjacentFace.Index]) { continue; }

                        isFaceModified[adjacentFace.Index] = true;

                        if ((i_Edge % 2) == 1) { i_NextFace_Weft.Add(adjacentFace.Index, i_Face); }
                        else { i_NextFace_Warp.Add(adjacentFace.Index, i_Face); }
                    }
                }

                foreach (int i_Face in i_CurrFace_Weft.Keys)
                {
                    // Move the face vertices.
                    List<HeEdge<Euc.Point>> edges = MoveVertices(deformedMesh, i_CurrFace_Weft[i_Face], i_Face, 1.0 / factor);

                    int nb_Edges = edges.Count;
                    for (int i_Edge = 1; i_Edge < nb_Edges; i_Edge++)
                    {
                        HeFace<Euc.Point> adjacentFace = edges[i_Edge].PairEdge.AdjacentFace;
                        if (adjacentFace == null || isFaceModified[adjacentFace.Index]) { continue; }

                        isFaceModified[adjacentFace.Index] = true;

                        if ((i_Edge % 2) == 1) { i_NextFace_Warp.Add(adjacentFace.Index, i_Face); }
                        else { i_NextFace_Weft.Add(adjacentFace.Index, i_Face); }
                    }
                }

                #endregion

            }

            #endregion

            #region Helpers

            /********** Primary Helpers **********/

            List<HeEdge<Euc.Point>> MoveVertices(HeMesh<Euc.Point> newMesh, int i_PrevFace, int i_Face, double factor)
            {
                // Get the sorted list of edges (starting at the edge common to i_Face and i_PrevFace)
                List<HeEdge<Euc.Point>> edges = SortedFaceEdges(newMesh, i_PrevFace, i_Face);

                // Get the sorted list of vertices (in correspondance with to the list of edges)
                List<HeVertex<Euc.Point>> vertices = new List<HeVertex<Euc.Point>>();
                foreach (HeEdge<Euc.Point> edge in edges) { vertices.Add(edge.StartVertex); }

                // If both vertices to move have already been modified, then leave the method
                if (isVertexModified[vertices[2].Index] & isVertexModified[vertices[3].Index]) { return edges; }

                // Get the normal vector of the face to build.
                Euc.Vector edge0 = new Euc.Vector(vertices[0].Position, vertices[1].Position);
                Euc.Vector normal = GetFaceNormal(newMesh, i_PrevFace, i_Face, -edge0, factor);

                // Modify Vertices
                if (!isVertexModified[vertices[2].Index])
                {
                    AssignVertexPosition(newMesh, normal, vertices[0].Index, vertices[1].Index, vertices[2].Index);
                }
                if (!isVertexModified[vertices[3].Index])
                {
                    AssignVertexPosition(newMesh, -normal, vertices[1].Index, vertices[0].Index, vertices[3].Index);
                }

                return edges;
            }

            /********** Secondary Helpers **********/

            // Returns the list of face edges which start at the edge common to i_Face and i_PrevFace
            List<HeEdge<Euc.Point>> SortedFaceEdges(HeMesh<Euc.Point> newMesh, int i_PrevFace, int i_Face)
            {
                List<HeEdge<Euc.Point>> sortedEdges = new List<HeEdge<Euc.Point>>();

                List<HeEdge<Euc.Point>> edges = newMesh.GetFace(i_Face).FaceEdges();
                int nb_Edges = edges.Count;

                int shift = 0;
                for (int i_Edge = 0; i_Edge < nb_Edges; i_Edge++)
                {
                    HeFace<Euc.Point> adjacentFace = edges[i_Edge].PairEdge.AdjacentFace;
                    if (!(adjacentFace is null) && adjacentFace.Index == i_PrevFace)
                    {
                        shift = i_Edge;
                        break;
                    }
                }

                for (int i_Edge = 0; i_Edge < nb_Edges; i_Edge++)
                {
                    HeEdge<Euc.Point> edge = edges[(i_Edge + shift) % nb_Edges];
                    sortedEdges.Add(edge);
                }

                return sortedEdges;
            }

            // Returns the normal vector of the face to build
            Euc.Vector GetFaceNormal(HeMesh<Euc.Point> newMesh, int i_PrevFace, int i_CurrFace, Euc.Vector commonEdge, double factor)
            {
                // Get the dihedral angle between movedFace and movingFace
                Euc.Vector old_NormalA = initialMesh.GetFace(i_PrevFace).Normal();
                Euc.Vector old_NormalB = initialMesh.GetFace(i_CurrFace).Normal();
                double dihedral = Euc.Vector.AngleBetween(old_NormalA, old_NormalB);

                List<HeEdge<Euc.Point>> oldEdges = SortedFaceEdges(initialMesh, i_PrevFace, i_CurrFace);
                Euc.Vector oldEdge0 = new Euc.Vector(oldEdges[0].StartVertex.Position, oldEdges[0].EndVertex.Position);
                Euc.Vector crossprod = Euc.Vector.CrossProduct(old_NormalA, old_NormalB);
                if (Euc.Vector.DotProduct(crossprod, oldEdge0) > 0) { dihedral = -dihedral; }

                dihedral = Math.Atan(factor * Math.Tan(dihedral));

                Euc.Vector new_NormalA = newMesh.GetFace(i_PrevFace).Normal();

                commonEdge.Unitize();
                new_NormalA.Rotate(new Euc.Point(0.0, 0.0, 0.0), commonEdge, /*factor * */ dihedral);

                new_NormalA.Unitize();

                return new_NormalA;
            }

            // Compute and assign the position of the VertexC in the newMesh
            void AssignVertexPosition(HeMesh<Euc.Point> newMesh, Euc.Vector normal, int i_VertexA, int i_VertexB, int i_VertexC)
            {
                // Get angle between the common edge and right edge (on the original mesh)
                Euc.Vector old_EdgeA = (Euc.Vector)(initialMesh.GetVertex(i_VertexA).Position - initialMesh.GetVertex(i_VertexB).Position);
                Euc.Vector old_EdgeB = (Euc.Vector)(initialMesh.GetVertex(i_VertexC).Position - initialMesh.GetVertex(i_VertexB).Position);

                double angle =  Euc.Vector.AngleBetween(old_EdgeA, old_EdgeB);

                // Get the right edge vector
                Euc.Vector edge = new Euc.Vector(newMesh.GetVertex(i_VertexB).Position, newMesh.GetVertex(i_VertexA).Position);
                
                edge.Rotate(new Euc.Point(0.0, 0.0, 0.0), normal, angle);
                edge.Unitize(); edge *= old_EdgeB.Length();

                // Assign vertex position
                newMesh.GetVertex(i_VertexC).Position = newMesh.GetVertex(i_VertexB).Position + edge;
                isVertexModified[i_VertexC] = true;
            }

            #endregion
        }
    }
}