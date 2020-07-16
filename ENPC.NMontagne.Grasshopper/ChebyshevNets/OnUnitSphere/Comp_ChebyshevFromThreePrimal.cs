using System;
using System.Collections.Generic;

using ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.ChebyshevNets;


namespace ENPC.NMontagne.Grasshopper.ChebyshevNet.OnUnitSphere
{
    /// <summary>
    /// A grasshopper component computing two patched Chebyshev net on the unity sphere from three primal conditions.
    /// </summary>
    public class Comp_ChebyshevFromThreePrimal : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_ChebyshevFromThreePrimal"/> class.
        /// </summary>
        public Comp_ChebyshevFromThreePrimal()
          : base("Chebyshev 3 Primal", "Cheb. 3P",
              "Computes two patched Chebyshev nets on the unit sphere from three primal conditions.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.ChebyshevNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Point(), "A Normals", "A", "The list of normals for the first Chebyshev net.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Point(), "B Normals", "B", "The list of normals common to the two Chebyshev nets.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Point(), "C Normals", "C", "The list of normals for the second Chebyshev net.", GH_K.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_HeMesh(), "HeMesh", "M", "The mesh representing the two patched Chebyshev nets.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            List<Point> a = new List<Point>();
            List<Point> b = new List<Point>();
            List<Point> c = new List<Point>(); 

            // Get Input

            if (!DA.GetDataList(0, a)) { return; }
            if (!DA.GetDataList(1, b)) { return; }
            if (!DA.GetDataList(2, c)) { return; }

            // Core of the component
            ChebyshevOnUnitSphere.Core_FromThreePrimal(a, b, c, out HeMesh<Point> mesh);

            // Set Output
            DA.SetData(0, mesh);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{8B14C890-E9F6-4FD1-8B19-9A094080B547}"); }
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
