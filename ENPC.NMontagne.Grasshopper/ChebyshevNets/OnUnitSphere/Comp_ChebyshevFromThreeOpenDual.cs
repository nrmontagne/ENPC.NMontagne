using System;
using System.Collections.Generic;

using ENPC.Geometry.Euclidean;
using ENPC.DataStructure.PolyhedralMesh.HalfedgeMesh;

using GH_K = Grasshopper.Kernel;

using ENPC.McNeel.Grasshopper.Parameters.Euclidean;

using ENPC.NMontagne.Core.CoreFunctions.ChebyshevNets;


namespace ENPC.NMontagne.Grasshopper.VossNets
{
    /// <summary>
    /// A grasshopper component computing two patched Chebyshev net on the unity sphere from three open dual conditions.
    /// </summary>
    public class Comp_ChebyshevFromThreeOpenDual : GH_K.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initializes a new instance of the <see cref="Comp_ChebyshevFromThreeOpenDual"/> class.
        /// </summary>
        public Comp_ChebyshevFromThreeOpenDual()
          : base("Chebyshev 2 Open Dual", "Cheb. 2OD",
              "Computes a Chebyshev net on the unit sphere from two open dual conditions.",
              LibrarySettings.Otter.Name, LibrarySettings.Otter.ChebyshevNets.Name)
        {
        }

        #endregion

        #region Override : GH_Component

        /// <inheritdoc cref="GH_K.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_K.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Param_Point(), "Main Diagonal A", "A0", "The list of normals of the main diagonal for the first Chebyshev.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Point(), "Minor Diagonal A", "A1", "The list of normals of a minor diagonal for the first Chebyshev.", GH_K.GH_ParamAccess.list);
            pManager.AddParameter(new Param_Point(), "Minor Diagonal B", "B1", "The mesh representing a the two patched Chebyshev nets.", GH_K.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_K.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_K.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Param_HeMesh(), "HeMesh", "M", "Mesh representing the Chebyshev net obtained from the input primal conditions.", GH_K.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_K.GH_Component.SolveInstance(GH_K.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_K.IGH_DataAccess DA)
        {
            // Declaration, initialization, and instanciation of input variables
            List<Point> a0 = new List<Point>();
            List<Point> a1 = new List<Point>();
            List<Point> b0 = new List<Point>();

            // Get Input

            if (!DA.GetDataList(0, a0)) { return; }
            if (!DA.GetDataList(1, a1)) { return; }
            if (!DA.GetDataList(2, b0)) { return; }

            // Core of the component
            ChebyshevOnUnitSphere.Core_FromThreeOpenDual(a0, a1, b0, out HeMesh<Point> mesh);

            // Set Output
            DA.SetData(0, mesh);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_K.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid
        {
            get { return new Guid("{E6BF4554-F899-4FCD-94C9-6E3E98FED494}"); }
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
