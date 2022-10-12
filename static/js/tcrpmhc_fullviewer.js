const DIR = tj.rundir_spath + "/" + jobid + "/";

var stage = new NGL.Stage("vwr");
stage.setParameters({backgroundColor: "white"});

window.addEventListener("resize", function(event) {
    stage.handleResize();
}, false);

stage.loadFile(modelfname).then(function(o) {
    o.autoView(1000);
    o.addRepresentation("cartoon", {
        name: "model",
        sele: "polymer",
        visible: true
    });
});

function loadCdrTmplt(file, name, resRange, chain) {
    var path = DIR + file;
    stage.loadFile(path).then(function(o) {
        o.autoView();
        o.addRepresentation("cartoon", {
            name: name,
            sele: resRange + ":" + chain,
            color: "lightgreen",
            visible: $("#"+name).prop("checked")
        });
    });
}

loadCdrTmplt(tj.acdr1_tmplt_id + "_Acdr1_tmplt.sup.pdb", "acdr1", "24-43", tj.acdr1_tmplt_id.slice(5,6));
loadCdrTmplt(tj.acdr2hv4_tmplt_id + "_Acdr2hv4_tmplt.sup.pdb", "acdr2", "56-91", tj.acdr2hv4_tmplt_id.slice(5,6));
loadCdrTmplt(tj.acdr3_tmplt_id + "_Acdr3_tmplt.sup.pdb", "acdr3", "107-139", tj.acdr3_tmplt_id.slice(5,6));
loadCdrTmplt(tj.bcdr1_tmplt_id + "_Bcdr1_tmplt.sup.pdb", "bcdr1", "24-43", tj.bcdr1_tmplt_id.slice(5,6));
loadCdrTmplt(tj.bcdr2hv4_tmplt_id + "_Bcdr2hv4_tmplt.sup.pdb", "bcdr2", "56-91", tj.bcdr2hv4_tmplt_id.slice(5,6));
loadCdrTmplt(tj.bcdr3_tmplt_id + "_Bcdr3_tmplt.sup.pdb", "bcdr3", "107-139", tj.bcdr3_tmplt_id.slice(5,6));

var path = DIR + tj.aori_tmplt_id + "_ori_tmplt.sup.pdb";
stage.loadFile(path).then(function(o) {
    o.autoView();
    o.addRepresentation("cartoon", {
        name: "orient",
        sele: "polymer",
        color: "pink",
        visible: $("#orient").prop("checked")
    });
});

$("#acdr1, #acdr2, #acdr3, #bcdr1, #bcdr2, #bcdr3, #orient").click(function() {
    var isChecked = $(this).prop("checked");
    var id = $(this).prop("id");
    stage.getRepresentationsByName(id).setVisibility(isChecked);
});

$("#centerbtn").click(function() {
  stage.autoView(1000);
});

$("#alignbtn").click(function() {
  var repr = stage.getRepresentationsByName("model").list[0];
  var axes = repr.repr.structure.getPrincipalAxes();
  stage.animationControls.rotate(axes.getRotationQuaternion(), 1500);
});
