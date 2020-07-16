using System;
using System.Collections.Generic;

using Euc = ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using Param_Euc = ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.VossNets;


namespace ENPC.NMontagne.Grasshopper.VossNets
{
    /// <summary>
    /// A grasshopper component computing a Voss net from a Gauss map using a target length for the edges.
    /// </summary>
    public class Comp_Voss_FromTwoPrimal : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_Voss_FromTwoPrimal"/> class.
        /// </summary>
        public Comp_Voss_FromTwoPrimal()
          : base("Voss 2 Primal", "V 2P.",
              "Computes a Voss net from two primal conditions.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.VossNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_Point(), "Normals U", "Nu", "The normals in u-direction.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Euc.Param_Point(), "Normals V", "Nv", "The normals in v-direction.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Euc.Param_Vector(), "Edges U", "Eu", "The edge vector in u-direction.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Euc.Param_Vector(), "Edges V", "Ev", "The edge vector in u-direction.", GH_K.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Gauss map", "G", "The Gauss map generated from the two primal conditions.", GH_K.GH_ParamAccess.item);
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Voss net", "V", "The Voss net generated from the two primal conditions.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            List<Euc.Point> normmals_U = new List<Euc.Point>();
            List<Euc.Point> normmals_V = new List<Euc.Point>();
            List<Euc.Vector> edges_U = new List<Euc.Vector>();
            List<Euc.Vector> edges_V = new List<Euc.Vector>();
           
            // Get Input

            if (!DA.GetDataList(0, normmals_U)) { return; }
            if (!DA.GetDataList(1, normmals_V)) { return; }
            if (!DA.GetDataList(2, edges_U)) { return; }
            if (!DA.GetDataList(3, edges_V)) { return; }

            // Core of the component
            DirectVoss.Core_FromTwoPrimal(normmals_U, normmals_V, edges_U, edges_V, out HeMesh<Euc.Point> gaussMap, out HeMesh<Euc.Point> vossNet);

            // Set Output
            DA.SetData(0, gaussMap);
            DA.SetData(0, vossNet);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{D61934B9-387C-4E71-A752-F997F054D847}"); }
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
