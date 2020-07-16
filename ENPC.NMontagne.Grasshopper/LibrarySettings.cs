using System;

namespace ENPC.NMontagne.Grasshopper
{
    internal static class LibrarySettings
    {
        #region Internal Variable

        /// <summary>
        /// Angular precision (in radians)
        /// </summary>
        internal const double _angularPrecision = Math.PI / 1e-5;

        /// <summary>
        /// Absolute linear
        /// </summary>
        internal const double _absolutePrecision = 1e-8;

        #endregion

        #region ENPC Grasshopper Library

        /// <summary>
        /// ENPC grashopper category.
        /// </summary>
        internal static class Otter
        {
            /// <summary>
            /// Name of the ENPC grasshopper library.
            /// </summary>
            public const string Name = "Otter";
            /// <summary>
            /// Nickname of the ENPC grasshopper library.
            /// </summary>
            public const string Nickname = "Otter";

            /************************* ENPC Subcategories *************************/

            /// <summary>
            /// ENPC grasshopper subcategory for Voss nets.
            /// </summary>
            internal static class Meshes
            {
                /// <summary>
                /// Name of the grasshopper subcategory.
                /// </summary>
                public const string Name = "Meshes";
                /// <summary>
                /// Nickname of the grasshopper subcategory.
                /// </summary>
                public const string Nickname = "Meshes";
            }

            /// <summary>
            /// ENPC grasshopper subcategory for Voss nets.
            /// </summary>
            internal static class VossNets
            {
                /// <summary>
                /// Name of the grasshopper subcategory.
                /// </summary>
                public const string Name = "Voss Nets";
                /// <summary>
                /// Nickname of the grasshopper subcategory.
                /// </summary>
                public const string Nickname = "Voss";
            }

            /// <summary>
            /// ENPC grasshopper subcategory for Chebyshev nets.
            /// </summary>
            internal static class ChebyshevNets
            {
                /// <summary>
                /// Name of the grasshopper subcategory.
                /// </summary>
                public const string Name = "Chebyshev Nets";
                /// <summary>
                /// Nickname of the grasshopper subcategory.
                /// </summary>
                public const string Nickname = "Cheb.";
            }
        }

        #endregion

    }
}
