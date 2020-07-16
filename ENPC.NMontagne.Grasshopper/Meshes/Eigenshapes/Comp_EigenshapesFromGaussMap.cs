using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using Param_Euc = ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.Meshes;
using Rhino.Geometry;


namespace ENPC.NMontagne.Grasshopper.Meshes
{
    /// <summary>
    /// A grasshopper component computing the eigenshapes of a net from its Gauss map.
    /// </summary>
    public class Comp_EigenshapesFromGaussMap : GH_K.GH_Component
    {
        #region Fields

        /// <summary>
        /// The mesh on which the EigenShapes were computed.
        /// </summary>
        public static HeMesh<Euc.Point> _mesh;
        /// <summary>
        /// The eigen modes of the stored mesh (Keys: mode index, Values: Directions for each vertices). 
        /// </summary>
        public static Dictionary<int, Euc.Vector[]> _modes;

        #endregion

        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_EigenshapesFromGaussMap"/> class.
        /// </summary>
        public Comp_EigenshapesFromGaussMap()
          : base("Eigenshapes G.", "Eig. G.",
              "Computes the eigenshapes of a net from its Gauss map.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.Meshes.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Gauss Map", "G", "The Gauss map of the Voss net to compute.", GH_K.GH_ParamAccess.item);
            pManager.AddIntegerParameter("Modes Indices", "I", "Defines the modes to apply.", GH_K.GH_ParamAccess.list);
            pManager.AddNumberParameter("Modes Amplitude", "A", "Defines the amplitude with which the modes are applied.", GH_K.GH_ParamAccess.list);
            pManager.AddBooleanParameter("Update Modes", "U", "Updates the eigenshapes according to the provided Gauss map.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Voss net", "V", "The Voss net generated from the input Gauss map.", GH_K.GH_ParamAccess.item);
            pManager.AddIntegerParameter("Number of Modes", "N", "Indicates the number of eigenshapes (or deformation modes) available.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            HeMesh<Euc.Point> gaussMap = new HeMesh<Euc.Point>();
            List<int> i_Modes = new List<int>();
            List<double> amplitudes = new List<double>();
            bool Update = false;

            // Get Input

            if (!DA.GetData(0, ref gaussMap)) { return; }
            if (!DA.GetDataList(1, i_Modes)) { return; }
            if (!DA.GetDataList(2, amplitudes)) { return; }
            if (!DA.GetData(3, ref Update)) { return; }

            // Core of the component
            if (Update)
            {
                Eigenshapes.Core_FromGaussMap(gaussMap, out _mesh, out _modes);
            }

            HeMesh<Euc.Point> otherMesh = new HeMesh<Euc.Point>();
            if (_mesh is null | _modes is null) { throw new NullReferenceException("The mesh or the modes were not initialized."); }
            else
            {
                otherMesh = (HeMesh<Euc.Point>)_mesh.Clone();

                int nb_Vertex = otherMesh.VertexCount;
                for (int index = 0; index < i_Modes.Count; index++)
                {
                    for (int i_Vertex = 0; i_Vertex < nb_Vertex; i_Vertex++)
                    {
                        otherMesh.GetVertex(i_Vertex).Position += (Euc.Point)(amplitudes[index] * _modes[i_Modes[index]][i_Vertex]);
                    }
                }
            }

            // Set Output
            DA.SetData(0, otherMesh);
            DA.SetData(1, _modes.Count);

        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{467BAC04-C381-4544-AA1E-0B4F50638BF8}"); }
        }

        /// <inheritdoc cref="GH_K.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <inheritdoc cref="GH_K.GH_DocumentObject.Exposure"/>
        public override GH_K.GH_Exposure Exposure
        {
            get { return GH_K.GH_Exposure.primary; }
        }

        #endregion
    }
}