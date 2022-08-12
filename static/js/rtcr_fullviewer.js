var stage = new NGL.Stage("vwr");
stage.setParameters({backgroundColor: "white"});

window.addEventListener("resize", function(event) {
    stage.handleResize();
}, false);

stage.loadFile(modelfname).then(function(o) {
    o.autoView();
    o.addRepresentation("cartoon", {
        name: "model",
        sele: "polymer",
        visible: true
    });
});

function loadCdrTmplt(file, name, resRange, chain) {
    var path = tj.rundir_spath + "/" + jobid + "/" + file;
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

loadCdrTmplt(tj.acdr1_tmplt_pdb + "_Acdr1_tmplt.pdb", "acdr1", "24-43", tj.acdr1_tmplt_pdb_chain);
loadCdrTmplt(tj.acdr2hv4_tmplt_pdb + "_Acdr2hv4_tmplt.pdb", "acdr2", "56-91", tj.acdr2hv4_tmplt_pdb_chain);
loadCdrTmplt(tj.acdr3_tmplt_pdb + "_Acdr3_tmplt.pdb", "acdr3", "107-139", tj.acdr3_tmplt_pdb_chain);
loadCdrTmplt(tj.bcdr1_tmplt_pdb + "_Bcdr1_tmplt.pdb", "bcdr1", "24-43", tj.bcdr1_tmplt_pdb_chain);
loadCdrTmplt(tj.bcdr2hv4_tmplt_pdb + "_Bcdr2hv4_tmplt.pdb", "bcdr2", "56-91", tj.bcdr2hv4_tmplt_pdb_chain);
loadCdrTmplt(tj.bcdr3_tmplt_pdb + "_Bcdr3_tmplt.pdb", "bcdr3", "107-139", tj.bcdr3_tmplt_pdb_chain);

function loadOriTmplt(file) {
    var path = tj.rundir_spath + "/" + jobid + "/" + file;
    stage.loadFile(path).then(function(o) {
        o.autoView();
        o.addRepresentation("cartoon", {
            name: "orient",
            sele: "polymer",
            color: "pink",
            visible: $("#orient").prop("checked")
        });
    });
}

loadOriTmplt(tj.ori_tmplt_Apdb + "_oriA_tmplt.pdb");
loadOriTmplt(tj.ori_tmplt_Bpdb + "_oriB_tmplt.pdb");

$("#acdr1, #acdr2, #acdr3, #bcdr1, #bcdr2, #bcdr3, #orient").click(function() {
    var isChecked = $(this).prop("checked");
    stage.getRepresentationsByName($(this).prop("id")).setVisibility(isChecked);
});

$("#centerbtn").click(function() {
  stage.autoView(1000);
});

$("#alignbtn").click(function() {
  var repr = stage.getRepresentationsByName("model").list[0];
  var axes = repr.repr.structure.getPrincipalAxes();
  stage.animationControls.rotate(axes.getRotationQuaternion(), 1500);
});
