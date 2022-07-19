// Create NGL Stage object
var stage = new NGL.Stage("nglviewport");
stage.setParameters({backgroundColor: "white"});

// handle window resizing
window.addEventListener("resize", function(event) {
    stage.handleResize();
}, false );

function loadStructure(input) {
  stage.removeAllComponents()
  return stage.loadFile(input).then(function(o) {
    o.autoView()
    o.addRepresentation($("#reprsdd :selected").val(), {
      sele: "polymer",
      name: "polymer",
      colorScheme: $("#colorsdd :selected").val()
    });
  });
}

function polymerRender(e) {
  stage.getRepresentationsByName("polymer").dispose()
  stage.eachComponent(function(o) {
    var selRepr = o.addRepresentation(e, {
      name: "polymer",
      sele: "polymer"
    });
    $("#colorsdd").val(selRepr.repr.colorScheme);
  });
}

$("#reprsdd").change(function() {
  var selRepr = $(this).val();
  polymerRender(selRepr);
});

$("#colorsdd").change(function() {
  var selColors = $(this).val();
  stage.getRepresentationsByName("polymer").setColor(selColors);
});

$("#spin-checkbox").click(function() {
  stage.toggleSpin();
});

$("#centerbtn").click(function() {
  stage.autoView(1000);
});

$("#alignbtn").click(function() {
  var repr = stage.getRepresentationsByName("polymer").list[0];
  var axes = repr.repr.structure.getPrincipalAxes();
  stage.animationControls.rotate(axes.getRotationQuaternion(), 1500);
});

loadStructure(modelfname);
