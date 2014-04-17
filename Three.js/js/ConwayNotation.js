var Conway = {

    buildPolyhedron: function (conwayNotation) {
        if (Conway.goodInput(conwayNotation)) {
            var polyhedron = Conway.generatePoly(conwayNotation);  // create poly.
            console.log(polyhedron);
            return  {
                "name": polyhedron.name,
                "category": ["?"],
                "vertex": polyhedron.xyz,
                "edge": Conway.buildEdges(polyhedron.face),
                "face": polyhedron.face
            };
        }
        return null;
    },

    buildEdges: function (faces) {
        var edgesDictionary = {};
        for (var i = 0; i < faces.length; i++) {
            var face = faces[i];
            for (var j = 0 ; j < face.length ; j++) {
                var v1 = face[j];
                var v2 = face[(j+1)%face.length];
                if (v2 < v1) {
                    v1 = v2;
                    v2 = face[j];
                }
                if (!edgesDictionary.hasOwnProperty(v1)) {
                    edgesDictionary[v1] = {};
                }
                edgesDictionary[v1][v2] = true;
            }
        }

        var edges = [];
        for (var firstEnd in edgesDictionary) {
            if (edgesDictionary.hasOwnProperty(firstEnd)) {
                for (var secondEnd in edgesDictionary[firstEnd]) {
                    if (edgesDictionary[firstEnd].hasOwnProperty(secondEnd)) {
                        edges.push([firstEnd, secondEnd]);
                    }
                }
            }
        }
        return edges;
    },

    goodInput: function (question) {                        // check input format
        if (question.search(/([^ktajsgebomdcrpTCOIDPAY0123456789])/) != -1) {
            console.error("Undefined character: " + RegExp.lastParen);
            return (false);
        }
        if (question.search(/^([ajsgebomdcrp]|([kt][0-9]*))*([TCOID]$|[PAY][0-9]+$)/) == -1) {
            console.error("Ill-formed polyhedron expression: " + question);
            return (false);
        }
        if (question.search(/([0-9]+)$/) != -1 && RegExp.lastParen < 3) {
            console.error("At least 3 sides are required, not " + RegExp.lastParen);
            return (false);
        }
        return (true);       // found no problems
    },

    state: function (x) {                         // put message in status bar
        window.status = x;
    },

    //------------------notation processing functions--------------------

    generatePoly: function (notation) {     // create polyhedron from notation
        var poly = new Conway.Polyhedron();       // each polyhedron, during construction
        var n = 0;                           // numeric argument

        Conway.state("Analyzing input...");
        var ops = Conway.getOps(notation);

        if (ops == Conway.globSavedPoly.name)   // if same as last time, use stored poly.
            return(Conway.globSavedPoly);
        if (ops == Conway.globSavedDual.name)   // if dual of last time, use stored dual.
            return(Conway.dual());
        if (ops.substr(-Conway.globSavedPoly.name.length) == Conway.globSavedPoly.name) {
            ops = ops.substr(0, ops.length - Conway.globSavedPoly.name.length);
            poly = Conway.globSavedPoly;         // extend previous poly
        }
        else if (ops.substr(-Conway.globSavedDual.name.length) == Conway.globSavedDual.name) {
            ops = ops.substr(0, ops.length - Conway.globSavedDual.name.length);
            poly = Conway.dual();                // extend previous dual
        }
        else {                           // start afresh
            if (ops.search(/([0-9]+)$/) != -1) {    // get number if present
                n = 1 * RegExp.lastParen;
                ops = ops.slice(0, -RegExp.lastParen.length);
            }
            Conway.state("Constructing seed...");
            if (ops.slice(-1) == "T") poly = Conway.tetrahedron();
            if (ops.slice(-1) == "O") poly = Conway.octahedron();
            if (ops.slice(-1) == "C") poly = Conway.cube();
            if (ops.slice(-1) == "I") poly = Conway.icosahedron();
            if (ops.slice(-1) == "D") poly = Conway.dodecahedron();
            if (ops.slice(-1) == "P") poly = Conway.prism(n);
            if (ops.slice(-1) == "A") poly = Conway.antiprism(n);
            if (ops.slice(-1) == "Y") poly = Conway.pyramid(n);
        }
        for (; ops != "";) {            // while loop
            n = 0;
            if (ops.search(/([0-9]+)$/) != -1) {    // get number if present
                n = 1 * RegExp.lastParen;
                ops = ops.slice(0, -RegExp.lastParen.length);
                if (n < 3)
                    console.warn("Of course you know that a value of " + n + " makes no sense, but I'll look anyway.");
            }
            if (ops.slice(-1) == "k") poly = Conway.kisN(poly, n);
            if (ops.slice(-1) == "a") poly = Conway.ambo(poly);
            if (ops.slice(-1) == "g") poly = Conway.gyro(poly);
            if (ops.slice(-1) == "p") poly = Conway.propellor(poly);
            if (ops.slice(-1) == "d") poly = Conway.dual();    // dual already computed
            if (ops.slice(-1) == "c") poly.xyz = Conway.canonicalXYZ(poly, 10);
            if (ops.slice(-1) == "r") poly = Conway.reflect(poly);
            ops = ops.slice(0, -1);   // remove last character
        }
//   poly.xyz = canonicalXYZ(poly, 5)     // refine final coords of poly and dual
        return (poly);
    },

    getOps: function (question) {    //  Convert notation into string of ops
        var ans = question;              // Global replacements in notation:
        ans = ans.replace(/P4$/g, "C");  // P4 --> C   (C is prism)
        ans = ans.replace(/A3$/g, "O");  // A3 --> O   (O is antiprism)
        ans = ans.replace(/Y3$/g, "T");  // Y3 --> T   (T is pyramid)
        ans = ans.replace(/e/g, "aa");   // e --> aa   (abbr. for explode)
        ans = ans.replace(/b/g, "ta");   // b --> ta   (abbr. for bevel)
        ans = ans.replace(/o/g, "jj");   // o --> jj   (abbr. for ortho)
        ans = ans.replace(/m/g, "kj");   // m --> kj   (abbr. for meta)
        ans = ans.replace(/t(\d*)/g, "dk$1d");  // t(n) --> dk(n)d  (dual operations)
        ans = ans.replace(/j/g, "dad");  // j --> dad  (dual operations)
        ans = ans.replace(/s/g, "dgd");  // s --> dgd  (dual operations)
        ans = ans.replace(/dd/g, "");    // dd --> null  (order 2)
        ans = ans.replace(/ad/g, "a");   // ad --> a   (a_ = ad_)
        ans = ans.replace(/gd/g, "g");   // gd --> g   (g_ = gd_)
        ans = ans.replace(/aY/g, "A");   // aY --> A   (interesting fact)
        ans = ans.replace(/dT/g, "T");   // dT --> T   (self-dual)
        ans = ans.replace(/gT/g, "D");   // gT --> D   (symm change)
        ans = ans.replace(/aT/g, "O");   // aT --> O   (symm change)
        ans = ans.replace(/dC/g, "O");   // dC --> O   (dual pair)
        ans = ans.replace(/dO/g, "C");   // dO --> C   (dual pair)
        ans = ans.replace(/dI/g, "D");   // dI --> D   (dual pair)
        ans = ans.replace(/dD/g, "I");   // dD --> I   (dual pair)
        ans = ans.replace(/aO/g, "aC");  // aO --> aC  (for uniqueness)
        ans = ans.replace(/aI/g, "aD");  // aI --> aD  (for uniqueness)
        ans = ans.replace(/gO/g, "gC");  // gO --> gC  (for uniqueness)
        ans = ans.replace(/gI/g, "gD");  // gI --> gD  (for uniqueness)
        console.info(question + " executed as " + ans);
        return (ans);
    },

    globSavedPoly: {"face": [], "xyz": [], "name": "null polyhedron"},  // global. poly from last canonicalization
    globSavedDual: {"face": [], "xyz": [], "name": "null polyhedron"},  // global. dual from last canonicalization

    Polyhedron: function () {       // constructor of initially null polyhedron
        this.face = [];   // array of faces.          face.length = # faces
        this.xyz = [];    // array of vertex coords.  xyz.length = # of vertices
        this.name = "null polyhedron"
    },


//------------------------polyhedra functions------------------

// Topology stored as set of "faces."  Each face is list of n 0-based vertex indices
// corresponding to one n-sided face.  Vertices listed clockwise as seen from outside.
// See cube below for example.  The term "flag" refers to a triple of a face index and
// two adjacent vertex indices, in clockwise order.

// Two global variables save the most recently constructed polyhedron and its dual
// in case the user asks for the same one again but with different output options.
// These are updated after canonicalization or taking the dual:

//-------------------------primative polyhedra-----------------
    data: function (poly) {         // informative string
        var nEdges = poly.face.length + poly.xyz.length - 2;   // E = V + F -2
        return("(" + poly.face.length + " faces, " + nEdges + " edges, " + poly.xyz.length + " vertices)");
    },

    tetrahedron: function () {
        var ans = new Conway.Polyhedron();
        ans.name = "T";
        ans.face = new Array(
            new Array(0, 1, 2), new Array(0, 2, 3), new Array(0, 3, 1), new Array(1, 3, 2));
        ans.xyz = new Array(
            new Array(1., 1., 1.), new Array(1., -1., -1.),
            new Array(-1., 1., -1.), new Array(-1., -1., 1.))
        return (ans)
    },

    octahedron: function () {
        var ans = new Conway.Polyhedron();
        ans.name = "O";
        ans.face = new Array(
            new Array(0, 1, 2), new Array(0, 2, 3), new Array(0, 3, 4), new Array(0, 4, 1),
            new Array(1, 4, 5), new Array(1, 5, 2), new Array(2, 5, 3), new Array(3, 5, 4));
        ans.xyz = new Array(
            new Array(0, 0, 1.414), new Array(1.414, 0, 0), new Array(0, 1.414, 0),
            new Array(-1.414, 0, 0), new Array(0, -1.414, 0), new Array(0, 0, -1.414))
        return (ans)
    },

    cube: function () {
        var ans = new Conway.Polyhedron();
        ans.name = "C";
        ans.face = new Array(
            new Array(3, 0, 1, 2), new Array(3, 4, 5, 0), new Array(0, 5, 6, 1),
            new Array(1, 6, 7, 2), new Array(2, 7, 4, 3), new Array(5, 4, 7, 6));
        ans.xyz = new Array(
            new Array(0.707, 0.707, 0.707), new Array(-0.707, 0.707, 0.707),
            new Array(-0.707, -0.707, 0.707), new Array(0.707, -0.707, 0.707),
            new Array(0.707, -0.707, -0.707), new Array(0.707, 0.707, -0.707),
            new Array(-0.707, 0.707, -0.707), new Array(-0.707, -0.707, -0.707))
        return (ans)
    },

    icosahedron: function () {
        var ans = new Conway.Polyhedron();
        ans.name = "I";
        ans.face = new Array(
            new Array(0, 1, 2), new Array(0, 2, 3), new Array(0, 3, 4), new Array(0, 4, 5),
            new Array(0, 5, 1), new Array(1, 5, 7), new Array(1, 7, 6), new Array(1, 6, 2),
            new Array(2, 6, 8), new Array(2, 8, 3), new Array(3, 8, 9), new Array(3, 9, 4),
            new Array(4, 9, 10), new Array(4, 10, 5), new Array(5, 10, 7), new Array(6, 7, 11),
            new Array(6, 11, 8), new Array(7, 10, 11), new Array(8, 11, 9), new Array(9, 11, 10));
        ans.xyz = new Array(
            new Array(0, 0, 1.176), new Array(1.051, 0, 0.526),
            new Array(0.324, 1., 0.525), new Array(-0.851, 0.618, 0.526),
            new Array(-0.851, -0.618, 0.526), new Array(0.325, -1., 0.526),
            new Array(0.851, 0.618, -0.526), new Array(0.851, -0.618, -0.526),
            new Array(-0.325, 1., -0.526), new Array(-1.051, 0, -0.526),
            new Array(-0.325, -1., -0.526), new Array(0, 0, -1.176))
        return (ans)
    },

    dodecahedron: function () {
        var ans = new Conway.Polyhedron();
        ans.name = "D";
        ans.face = new Array(
            new Array(0, 1, 4, 7, 2), new Array(0, 2, 6, 9, 3), new Array(0, 3, 8, 5, 1),
            new Array(1, 5, 11, 10, 4), new Array(2, 7, 13, 12, 6), new Array(3, 9, 15, 14, 8),
            new Array(4, 10, 16, 13, 7), new Array(5, 8, 14, 17, 11), new Array(6, 12, 18, 15, 9),
            new Array(10, 11, 17, 19, 16), new Array(12, 13, 16, 19, 18), new Array(14, 15, 18, 19, 17));
        ans.xyz = new Array(
            new Array(0, 0, 1.07047), new Array(0.713644, 0, 0.797878),
            new Array(-0.356822, 0.618, 0.797878), new Array(-0.356822, -0.618, 0.797878),
            new Array(0.797878, 0.618034, 0.356822), new Array(0.797878, -0.618, 0.356822),
            new Array(-0.934172, 0.381966, 0.356822), new Array(0.136294, 1., 0.356822),
            new Array(0.136294, -1., 0.356822), new Array(-0.934172, -0.381966, 0.356822),
            new Array(0.934172, 0.381966, -0.356822), new Array(0.934172, -0.381966, -0.356822),
            new Array(-0.797878, 0.618, -0.356822), new Array(-0.136294, 1., -0.356822),
            new Array(-0.136294, -1., -0.356822), new Array(-0.797878, -0.618034, -0.356822),
            new Array(0.356822, 0.618, -0.797878), new Array(0.356822, -0.618, -0.797878),
            new Array(-0.713644, 0, -0.797878), new Array(0, 0, -1.07047))
        return (ans)
    },

    prism: function (n) {
        var theta = 6.283185 / n;        // pie angle
        var h = Math.sin(theta / 2);     // half-edge
        var ans = new Conway.Polyhedron();
        ans.name = "P" + n;

        for (var i1 = 0; i1 < n; i1++)     // vertex #'s 0...n-1 around one face
            ans.xyz[ans.xyz.length] = new Array(Math.cos(i1 * theta), Math.sin(i1 * theta), h);
        for (var i2 = 0; i2 < n; i2++)     // vertex #'s n...2n-1 around other
            ans.xyz[ans.xyz.length] = new Array(Math.cos(i2 * theta), Math.sin(i2 * theta), -h);

        ans.face[ans.face.length] = Conway.sequence(n - 1, 0);      // top
        ans.face[ans.face.length] = Conway.sequence(n, 2 * n - 1);    // bottom
        for (var i3 = 0; i3 < n; i3++)                            // n square sides:
            ans.face[ans.face.length] = new Array(i3, (i3 + 1) % n, (i3 + 1) % n + n, i3 + n);

        ans.xyz = Conway.adjustXYZ(ans, 1);
        return (ans);
    },

    antiprism: function (n) {
        var theta = 6.283185 / n;        // pie angle
        var h = Math.sqrt(1 - 4 / (4 + 2 * Math.cos(theta / 2) - 2 * Math.cos(theta))); // half-height
        var r = Math.sqrt(1 - h * h);      // radius of face circle
        var f = Math.sqrt(h * h + Math.pow(r * Math.cos(theta / 2), 2));
        r = r / f;  // correction so edge midpoints (not vertices) on unit sphere
        h = h / f;
        var ans = new Conway.Polyhedron();
        ans.name = "A" + n;

        for (var i1 = 0; i1 < n; i1++)     // vertex #'s 0...n-1 around one face
            ans.xyz[ans.xyz.length] = new Array(r * Math.cos(i1 * theta), r * Math.sin(i1 * theta), h);
        for (var i2 = 0; i2 < n; i2++)     // vertex #'s n...2n-1 around other
            ans.xyz[ans.xyz.length] = new Array(r * Math.cos((i2 + 0.5) * theta), r * Math.sin((i2 + 0.5) * theta), -h);

        ans.face[ans.face.length] = Conway.sequence(n - 1, 0);      // top
        ans.face[ans.face.length] = Conway.sequence(n, 2 * n - 1);    // bottom
        for (var i3 = 0; i3 < n; i3++) {                          // 2n triangular sides:
            ans.face[ans.face.length] = new Array(i3, (i3 + 1) % n, i3 + n);
            ans.face[ans.face.length] = new Array(i3, i3 + n, ((n + i3 - 1) % n + n));
        }
        ans.xyz = Conway.adjustXYZ(ans, 1);
        return (ans);
    },

    pyramid: function (n) {
        var theta = 6.283185 / n;        // pie angle
        var ans = new Conway.Polyhedron();
        ans.name = "Y" + n;

        for (var i1 = 0; i1 < n; i1++)     // vertex #'s 0...n-1 around base
            ans.xyz[ans.xyz.length] = new Array(Math.cos(i1 * theta), Math.sin(i1 * theta), .2);
        ans.xyz[ans.xyz.length] = new Array(0, 0, -2);    // apex

        ans.face[ans.face.length] = Conway.sequence(n - 1, 0);      // base
        for (var i2 = 0; i2 < n; i2++)                            // n triangular sides:
            ans.face[ans.face.length] = new Array(i2, (i2 + 1) % n, n);

        ans.xyz = Conway.canonicalXYZ(ans, 3);
        return (ans);
    },


    //----------------polyhedron operators---------------------------
    // Process: call newPoly() to clear tables
    //          for each vertex of new polyhedron:
    //              call newV(Vname, xyz) with symbolic name and approx location
    //          for each flag of new polyhedron:
    //              call newFlag(Fname, Vname1, Vname2)  with symbolic names
    //          call flags2poly()  to assemble flags into polyhedron structure
    //          canonicalize vertex locations
    //          set name as appropriate

    kisN: function (poly, n) {     // only kis n-sided faces, but n==0 means kiss all.
        Conway.state("Taking kis of " + (n == 0 ? "" : n + "-sided faces of ") + poly.name + "...");
        Conway.newPoly();
        for (var i1 = 0; i1 < poly.xyz.length; i1++)
            Conway.newV("v" + i1, poly.xyz[i1]);              // each old vertex is a new vertex
        var centers = Conway.faceCenters(poly);           // new vertices in centers of n-sided face
        var foundAny = false;                      // alert if don't find any
        for (var i2 = 0; i2 < poly.face.length; i2++) {
            var v1 = "v" + poly.face[i2][poly.face[i2].length - 1];  // previous vertex
            for (var j = 0; j < poly.face[i2].length; j++) {
                var v2 = "v" + poly.face[i2][j];                  // this vertex
                if (poly.face[i2].length == n || n == 0) {    // kiss the n's, or all
                    foundAny = true;                // flag that we found some
                    Conway.newV("f" + i2, centers[i2]);        // new vertex in face center
                    var fname = i2 + v1;
                    Conway.newFlag(fname, v1, v2);         // three new flags, if n-sided
                    Conway.newFlag(fname, v2, "f" + i2);
                    Conway.newFlag(fname, "f" + i2, v1);
                }
                else
                    Conway.newFlag(i2, v1, v2);             // same old flag, if non-n
                v1 = v2;                           // current becomes previous
            }
        }
        if (!foundAny)
            console.error("No " + n + "-fold components were found.");
        var ans = Conway.flags2poly();
        ans.name = "k" + (n == 0 ? "" : n) + poly.name;
        ans.xyz = Conway.adjustXYZ(ans, 3);               // adjust and
//   ans.xyz = canonicalXYZ(ans, 3);            // canonicalize lightly
        return (ans);
    },

    ambo: function (poly) {                      // compute ambo of argument
        Conway.state("Taking ambo of " + poly.name + "...");
        Conway.newPoly();
        for (var i = 0; i < poly.face.length; i++) {
            var v1 = poly.face[i][poly.face[i].length - 2];  // preprevious vertex
            var v2 = poly.face[i][poly.face[i].length - 1];  // previous vertex
            for (var j = 0; j < poly.face[i].length; j++) {
                var v3 = poly.face[i][j];        // this vertex
                if (v1 < v2)                     // new vertices at edge midpoints
                    Conway.newV(Conway.midName(v1, v2), Conway.midpoint(poly.xyz[v1], poly.xyz[v2]));
                Conway.newFlag("f" + i, Conway.midName(v1, v2), Conway.midName(v2, v3));     // two new flags
                Conway.newFlag("v" + v2, Conway.midName(v2, v3), Conway.midName(v1, v2));
                v1 = v2;                         // shift over one
                v2 = v3;
            }
        }
        var ans = Conway.flags2poly();
        ans.name = "a" + poly.name;
        ans.xyz = Conway.adjustXYZ(ans, 2);             // canonicalize lightly
        return (ans);
    },

    midName: function (v1, v2) {              // unique symbolic name, e.g. "1_2"
        if (v1 < v2)
            return (v1 + "_" + v2);
        else
            return (v2 + "_" + v1);
    },

    gyro: function (poly) {                      // compute gyro of argument
        Conway.state("Taking gyro of " + poly.name + "...");
        Conway.newPoly();
        for (var i1 = 0; i1 < poly.xyz.length; i1++)
            Conway.newV("v" + i1, Conway.unit(poly.xyz[i1]));           // each old vertex is a new vertex
        var centers = Conway.faceCenters(poly);              // new vertices in center of each face
        for (var i2 = 0; i2 < poly.face.length; i2++)
            Conway.newV("f" + i2, Conway.unit(centers[i2]));
        for (var i3 = 0; i3 < poly.face.length; i3++) {
            var v1 = poly.face[i3][poly.face[i3].length - 2];  // preprevious vertex
            var v2 = poly.face[i3][poly.face[i3].length - 1];  // previous vertex
            for (var j = 0; j < poly.face[i3].length; j++) {
                var v3 = poly.face[i3][j];                  // this vertex
                Conway.newV(v1 + "~" + v2, Conway.oneThird(poly.xyz[v1], poly.xyz[v2]));  // new v in face
                var fname = i3 + "f" + v1;
                Conway.newFlag(fname, "f" + i3, v1 + "~" + v2);          // five new flags
                Conway.newFlag(fname, v1 + "~" + v2, v2 + "~" + v1);
                Conway.newFlag(fname, v2 + "~" + v1, "v" + v2);
                Conway.newFlag(fname, "v" + v2, v2 + "~" + v3);
                Conway.newFlag(fname, v2 + "~" + v3, "f" + i3);
                v1 = v2;                                   // shift over one
                v2 = v3;
            }
        }
        var ans = Conway.flags2poly();
        ans.name = "g" + poly.name;
        ans.xyz = Conway.adjustXYZ(ans, 3);                       // canonicalize lightly
        return (ans);
    },

    propellor: function (poly) {                             // compute propellor of argument
        Conway.state("Taking propellor of " + poly.name + "...");
        Conway.newPoly();
        for (var i1 = 0; i1 < poly.xyz.length; i1++)
            Conway.newV("v" + i1, Conway.unit(poly.xyz[i1]));           // each old vertex is a new vertex
        for (var i2 = 0; i2 < poly.face.length; i2++) {
            var v1 = poly.face[i2][poly.face[i2].length - 2];  // preprevious vertex
            var v2 = poly.face[i2][poly.face[i2].length - 1];  // previous vertex
            for (var j = 0; j < poly.face[i2].length; j++) {
                var v3 = poly.face[i2][j];                  // this vertex
                Conway.newV(v1 + "~" + v2, Conway.oneThird(poly.xyz[v1], poly.xyz[v2]));  // new v in face
                var fname = i2 + "f" + v2;
                Conway.newFlag("v" + i2, v1 + "~" + v2, v2 + "~" + v3);      // five new flags
                Conway.newFlag(fname, v1 + "~" + v2, v2 + "~" + v1);
                Conway.newFlag(fname, v2 + "~" + v1, "v" + v2);
                Conway.newFlag(fname, "v" + v2, v2 + "~" + v3);
                Conway.newFlag(fname, v2 + "~" + v3, v1 + "~" + v2);
                v1 = v2;                                   // shift over one
                v2 = v3;
            }
        }
        var ans = Conway.flags2poly();
        ans.name = "p" + poly.name;
        ans.xyz = Conway.adjustXYZ(ans, 3);                       // canonicalize lightly
        return (ans);
    },

    reflect: function (poly) {                              // compute reflection through origin
        Conway.state("Taking reflection of " + poly.name + "...");
        for (var i1 = 0; i1 < poly.xyz.length; i1++)
            poly.xyz[i1] = Conway.mult(-1, poly.xyz[i1]);           // reflect each point
        for (var i2 = 0; i2 < poly.face.length; i2++)
            poly.face[i2] = poly.face[i2].reverse();         // repair clockwise-ness
        poly.name = "r" + poly.name;
        poly.xyz = Conway.adjustXYZ(poly, 1);                     // build dual
        return (poly);
    },

    //--------------------------------Dual------------------------------------------
    // the makeDual function computes the dual's topology, needed for canonicalization,
// where xyz's are determined.  It is then saved in a global variable globSavedDual.
    // when the d operator is executed, d just returns the saved value.

    dual: function () {        // for d operator, just swap poly with saved dual
        var ans = Conway.globSavedDual;
        Conway.globSavedDual = Conway.globSavedPoly;
        Conway.globSavedPoly = ans;
        return (ans);
    },

    makeDual: function (poly) {   // compute dual of argument, matching V and F indices
        Conway.state("Taking dual of " + poly.name + "...");
        Conway.newPoly();
        var face = [];            // make table of face as fn of edge
        for (var i1 = 0; i1 < poly.xyz.length; i1++)
            face[i1] = {};    // create empty associative table
        for (var i2 = 0; i2 < poly.face.length; i2++) {
            var vA1 = poly.face[i2][poly.face[i2].length - 1];  // previous vertex
            for (j = 0; j < poly.face[i2].length; j++) {
                var vA2 = poly.face[i2][j];                  // this vertex
                face[vA1]["v" + vA2] = i2;    // fill it.  2nd index is associative
                vA1 = vA2;                                   // current becomes previous
            }
        }
        for (var i3 = 0; i3 < poly.face.length; i3++)         // create d's v's per p's f's
            Conway.newV(i3, []);                      // only topology needed for canonicalize
        for (var i4 = 0; i4 < poly.face.length; i4++) {       // one new flag for each old one
            var vB1 = poly.face[i4][poly.face[i4].length - 1];  // previous vertex
            for (j = 0; j < poly.face[i4].length; j++) {
                var vB2 = poly.face[i4][j];                  // this vertex
                Conway.newFlag(vB1, face[vB2]["v" + vB1], i4);        // look up face across edge
                vB1 = vB2;                                   // current becomes previous
            }
        }
        var ans = Conway.flags2poly();      // this gives one indexing of answer
        var sortF = [];     // but f's of dual are randomly ordered, so sort
        for (var i = 0; i < ans.face.length; i++) {
            var j = Conway.intersect(poly.face[ans.face[i][0]], poly.face[ans.face[i][1]], poly.face[ans.face[i][2]]);
            sortF[j] = ans.face[i];  // p's v for d's f is common to three of p's f's
        }
        ans.face = sortF;            // replace with the sorted list of faces
        if (poly.name.substr(0, 1) != "d")
            ans.name = "d" + poly.name;        // dual name is same with "d" added...
        else
            ans.name = poly.name.substr(1);    // ...or removed
        return (ans);
    },

//-------------------Canonicalization Algorithm--------------------------
// True canonicalization rather slow.  Using center of gravity of vertices for each
// face gives a quick "adjustment" which planarizes faces at least.

    canonicalXYZ: function (poly, nIterations) {      // compute new vertex coords.
        var dpoly = Conway.makeDual(poly)     // v's of dual are in order or arg's f's
        Conway.state("Canonicalizing " + poly.name + "...");
        for (var count = 0; count < nIterations; count++) {    // iteration:
            dpoly.xyz = Conway.reciprocalN(poly);                 // reciprocate face normals
            poly.xyz = Conway.reciprocalN(dpoly);                 // reciprocate face normals
        }
        Conway.globSavedPoly = poly;      // save poly in global variable
        Conway.globSavedDual = dpoly;     // save dual in global variable
        return (poly.xyz);
    },

    reciprocalN: function (poly) {    // make array of vertices reciprocal to given planes
        var ans = [];
        for (var i = 0; i < poly.face.length; i++) {    // for each face:
            var centroid = Conway.vecZero();           // running sum of vertex coords
            var normal = Conway.vecZero();             // running sum of normal vectors
            var avgEdgeDist = 0.;               // running sum for avg edge distance
            var v1 = poly.face[i][poly.face[i].length - 2];  // preprevious vertex
            var v2 = poly.face[i][poly.face[i].length - 1];  // previous vertex
            for (var j = 0; j < poly.face[i].length; j++) {
                var v3 = poly.face[i][j];                  // this vertex
                centroid = Conway.add(centroid, poly.xyz[v3]);
                normal = Conway.add(normal, Conway.orthogonal(poly.xyz[v1], poly.xyz[v2], poly.xyz[v3]));
                avgEdgeDist = avgEdgeDist + Conway.edgeDist(poly.xyz[v1], poly.xyz[v2]);
                v1 = v2;                                   // shift over one
                v2 = v3;
            }
            centroid = Conway.mult(1 / poly.face[i].length, centroid);
            normal = Conway.unit(normal);
            avgEdgeDist = avgEdgeDist / poly.face[i].length;
            ans[i] = Conway.reciprocal(Conway.mult(Conway.dot(centroid, normal), normal));  // based on face
            ans[i] = Conway.mult((1 + avgEdgeDist) / 2, ans[i]);                  // edge correction
        }
        return (ans);
    },

    adjustXYZ: function (poly, nIterations) {           // quick planarization
        var dpoly = Conway.makeDual(poly)     // v's of dual are in order or arg's f's
        Conway.state("Planarizing " + poly.name + "...");
        for (var count = 0; count < nIterations; count++) {    // iteration:
            dpoly.xyz = Conway.reciprocalC(poly);             // reciprocate face centers
            poly.xyz = Conway.reciprocalC(dpoly);             // reciprocate face centers
        }
        Conway.globSavedPoly = poly;      // save poly in global variable
        Conway.globSavedDual = dpoly;     // save dual in global variable
        return (poly.xyz);
    },

    reciprocalC: function (poly) {           // return array of reciprocals of face centers
        var center = Conway.faceCenters(poly);
        for (var i = 0; i < poly.face.length; i++) {
            var m2 = center[i][0] * center[i][0] + center[i][1] * center[i][1] + center[i][2] * center[i][2];
            center[i][0] = center[i][0] / m2;   // divide each coord by magnitude squared
            center[i][1] = center[i][1] / m2;
            center[i][2] = center[i][2] / m2;
        }
        return(center);
    },

    faceCenters: function (poly) {              // return array of "face centers"
        var ans = [];
        for (var i = 0; i < poly.face.length; i++) {
            ans[i] = Conway.vecZero();                      // running sum
            for (var j = 0; j < poly.face[i].length; j++)    // just average vertex coords:
                ans[i] = Conway.add(ans[i], poly.xyz[poly.face[i][j]]);  // sum and...
            ans[i] = Conway.mult(1. / poly.face[i].length, ans[i]);        // ...divide by n
        }
        return (ans);
    },

    //----------------polyhedron assembly from flags-------------------
    // 4 global objects used, since javascript won't pass by reference
// property lists used as associative arrays of symbolic names

    gFlag: {},   // gFlag[face][vertex]=next vertex of flag; symbolic triples
    gXYZ: {},   // XYZ coordinates
    gVert: {},   // [symbolic names] holds vertex index
    gFace: {},   // list of symbolic names for faces

    newPoly: function () {     // clear global vars in preparation for new construction
        Conway.gFlag = {};
        Conway.gXYZ = {};
        Conway.gVert = {};
        Conway.gFace = {};
    },

    newFlag: function (face, v1, v2) {    // add flag and face to list
        if (Conway.gFlag[face] == null)
            Conway.gFlag[face] = {};  // create entry for face if needed
        Conway.gFlag[face][v1] = v2;            // create next-vertex entry
    },

    newV: function (name, xyz) {      // add vertex, if new, to lists
        if (Conway.gVert[name] == null) {
            Conway.gVert[name] = 0;         // dummy value for now
            Conway.gXYZ[name] = xyz;
        }
    },

    flags2poly: function () {     // arrange symbolic flags into polyhedron format
        var poly = new Conway.Polyhedron();
        var ctr = 0;                     // first number the vertices
        for (var i1 in Conway.gVert) {
            if (Conway.gVert.hasOwnProperty(i1)) {
                poly.xyz[ctr] = Conway.gXYZ[i1];   // and store in array
                Conway.gVert[i1] = ctr;
                ctr++;
            }
        }
        ctr = 0;                    // now number the faces
        for (var i2 in Conway.gFlag) {      // for each face
            if (Conway.gFlag.hasOwnProperty(i2)) {
                poly.face[ctr] = [];
                var v0;                  // any vertex as starting point
                for (var j in Conway.gFlag[i2]) {
                    if (Conway.gFlag[i2].hasOwnProperty(j)) {
                        v0 = Conway.gFlag[i2][j];
                        break;               // need just one.
                    }
                }
                var v = v0;              // v moves around face
                do {
                    poly.face[ctr][poly.face[ctr].length] = Conway.gVert[v];   // record index
                    v = Conway.gFlag[i2][v];                                // go to next vertex
                } while (v != v0);                              // until back to start
                ctr++;
            }
        }
        Conway.newPoly();                  // release memory
        poly.name = "unknown polyhedron"
        return (poly);
    },

//-----------------------math functions--------------------------

    vecZero: function () {
        var ans = [];
        ans[0] = 0.;
        ans[1] = 0.;
        ans[2] = 0.;
        return (ans);
    },

    mult: function (c, vec) {       // c times 3-vector
        var ans = [];
        ans[0] = c * vec[0];
        ans[1] = c * vec[1];
        ans[2] = c * vec[2];
        return (ans);
    },

    add: function (vec1, vec2) {    // sum two 3-vectors
        var ans = [];
        ans[0] = vec1[0] + vec2[0];
        ans[1] = vec1[1] + vec2[1];
        ans[2] = vec1[2] + vec2[2];
        return (ans);
    },

    sub: function (vec1, vec2) {    // subtract two 3-vectors
        var ans = [];
        ans[0] = vec1[0] - vec2[0];
        ans[1] = vec1[1] - vec2[1];
        ans[2] = vec1[2] - vec2[2];
        return (ans);
    },

    dot: function (vec1, vec2) {    // dot product two 3-vectors
        return (vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]);
    },

    midpoint: function (vec1, vec2) {    // mean of two 3-vectors
        var ans = [];
        ans[0] = 0.5 * (vec1[0] + vec2[0]);
        ans[1] = 0.5 * (vec1[1] + vec2[1]);
        ans[2] = 0.5 * (vec1[2] + vec2[2]);
        return (ans);
    },

    oneThird: function (vec1, vec2) {    // approx. (2/3)v1 + (1/3)v2   (assumes 3-vector)
        var ans = [];
        ans[0] = 0.7 * vec1[0] + 0.3 * vec2[0];
        ans[1] = 0.7 * vec1[1] + 0.3 * vec2[1];
        ans[2] = 0.7 * vec1[2] + 0.3 * vec2[2];
        return (ans);
    },

    reciprocal: function (vec) {    // reflect 3-vector in unit sphere
        var factor = 1. / Conway.mag2(vec);
        var ans = [];
        ans[0] = factor * vec[0];
        ans[1] = factor * vec[1];
        ans[2] = factor * vec[2];
        return (ans);
    },

    unit: function (vec) {          // normalize 3-vector to unit magnitude
        var size = Conway.mag2(vec);
        if (size == 0.) {          // remove this test someday...
            console.error("Mag(zero) -- probable bug.");
            return (vec);
        }
        var c = 1. / Math.sqrt(size);
        var ans = [];
        ans[0] = c * vec[0];
        ans[1] = c * vec[1];
        ans[2] = c * vec[2];
        return (ans);
    },

    mag2: function (vec) {          // magnitude squared of 3-vector
        return (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    },

    tangentPoint: function (v1, v2) {    // point where line v1...v2 tangent to an origin sphere
        var d = Conway.sub(v2, v1);             // difference vector
        return(Conway.sub(v1, Conway.mult(Conway.dot(d, v1) / Conway.mag2(d), d)));
    },

    edgeDist: function (v1, v2) {         // distance of line v1...v2 to origin
        return(Math.sqrt(Conway.mag2(Conway.tangentPoint(v1, v2))));
    },

    orthogonal: function (v3, v2, v1) {   // find unit vector orthog to plane of 3 pts
        var d1 = Conway.sub(v2, v1);            // adjacent edge vectors
        var d2 = Conway.sub(v3, v2);
        var ans = [];
        ans[0] = d1[1] * d2[2] - d1[2] * d2[1];    // cross product
        ans[1] = d1[2] * d2[0] - d1[0] * d2[2];
        ans[2] = d1[0] * d2[1] - d1[1] * d2[0];
        return (ans)
    },

    intersect: function (set1, set2, set3) {  // find element common to 3 sets
        for (var i = 0; i < set1.length; i++)    // by brute force search
            for (var j = 0; j < set2.length; j++)
                if (set1[i] == set2[j])
                    for (var k = 0; k < set3.length; k++)
                        if (set1[i] == set3[k])
                            return (set1[i]);
        console.error("program bug in intersect()");
        return (null);
    },

    sequence: function (start, stop) {    // make list of integers, inclusive
        var ans = [];
        if (start <= stop)
            for (var i1 = start; i1 <= stop; i1++)
                ans[ans.length] = i1;
        else
            for (var i2 = start; i2 >= stop; i2--)
                ans[ans.length] = i2;
        return (ans);
    }

}
