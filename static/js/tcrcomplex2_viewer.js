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
    o.addRepresentation("ball+stick", {
      name: "ligand",
      visible: $("#ligand-checkbox").checked,
      sele: "not ( polymer or water or ion )"
    });
    o.addRepresentation("spacefill", {
      name: "waterIon",
      visible: $("#water-ion-checkbox").checked,
      sele: "water or ion",
      scale: 0.25
    });
  });
}

function polymerRender(e) {
  stage.getRepresentationsByName("polymer").dispose()
  stage.eachComponent(function(o) {
    var selRepr = o.addRepresentation(e, {
      sele: "polymer",
      name: "polymer",
    });
    $("#colorsdd").val(selRepr.repr.colorScheme);
  });
}

$("#reprsdd").change(function() {
  var selRepr = $(this).val();
  polymerRender(selRepr);
});

$("#ligand-checkbox").click(function() {
  var isChecked = $(this).checked;
  stage.getRepresentationsByName("ligand").setVisibility(isChecked);
});

$("#water-ion-checkbox").click(function() {
  var isChecked = $(this).checked;
  stage.getRepresentationsByName("waterIon").setVisibility(isChecked);
});

$("#colorsdd").change(function() {
  var selColors = $(this).val();
  stage.getRepresentationsByName("polymer").setColor(selColors);
});

$("#spin-checkbox").click(function() {
  stage.toggleSpin();
});

$("#centerbtn").click(function(e) {
  e.preventDefault();
  stage.autoView(1000);
});

loadStructure(modelfname);
