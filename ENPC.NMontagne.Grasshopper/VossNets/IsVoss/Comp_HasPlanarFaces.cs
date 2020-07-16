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
    /// A grasshopper component evaluating whether the net has planar faces.
    /// </summary>
    public class Comp_HasPlanarFaces : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_HasPlanarFaces"/> class.
        /// </summary>
        public Comp_HasPlanarFaces()
          : base("Is Planar", "Is Plan.",
              "Evaluating whether the net has planar faces.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.VossNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_HeMesh(), "Net", "M", "The net to evaluate.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_Euc.Param_Point(), "Valid Faces", "T", "The barycenter of the faces for which the planarity criteria is valid.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Euc.Param_Point(), "Invalid Faces", "F", "The barycenter of the faces for which the planarity criteria is not valid.", GH_K.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            HeMesh<Euc.Point> net = new HeMesh<Euc.Point>();

            // Get Input

            if (!DA.GetData(0, ref net)) { return; }

            // Core of the component
            IsVoss.Core_HasPlanarFaces(net, out List<Euc.Point> areTrue, out List<Euc.Point> areFalse);

            // Set Output
            DA.SetDataList(0, areTrue);
            DA.SetDataList(1, areFalse);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{21C62099-C8EB-4E32-A9F1-72664BFD6AC8}"); }
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
