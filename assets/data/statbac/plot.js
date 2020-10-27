var globals = {};

var split_by_params = {
    title: "Admission au bac (% d'inscrits)",
    description: "% d'admis au baccalauréat depuis 1969",
    full_width: true,
    height: 300,
    bottom: 65,
    left: 50,
    right: 30,
    xax_count: 4,
    missing_is_hidden: true,
    missing_text : 'Inscription|Admission au bac depuis 1969',
    target: '#bacperyear',
    x_accessor: 'years',
    yax_units: '%',
    show_tooltips: true,
    y_accessor: 'adm',
    x_label: 'années',
    mouseover: function(d, i) {
            // custom format the rollover text, show days
            d3.select('#bacperyear svg .mg-active-datapoint')
                .html(d.cand !=null ? '- y: '+ d.years + ' | ' + d.cand +  ' ( ' + d.adm + '% )': '- Année blanche -');
        }


};

var resstagne = {
    title: "Prédiction du % d'admission au bac",
    description: "Prédiction du pourcentage d'admission au baccalauréat jusqu'en 2020, en utilisant la regréssion linéaire avec RANSAC.",
    full_width: true,
    height: 350,
    bottom: 65,
    left: 50,
    right: 60,
    area: false,
    missing_is_hidden: true,
    yax_units: '%',
    target: '#resstagne',
    x_accessor: 'years',
    show_tooltips: true,
    y_accessor: ['adm', 'adm_pred'],
    x_label: 'années',
    aggregate_rollover: true,
    baselines: [{value: 33.33, label: '33.3%, la limite'}],
    legend : ['réel', 'prédict. (Lin. regr.)']
};


var primary_completion = {
    title: "Taux d'achèvement du primaire",
    description: "Taux d'achèvement de l'école primaire au Bénin par année",
    full_width: true,
    height: 350,
    bottom: 65,
    left: 50,
    right: 30,
    area: false,
    xax_count: 4,
    missing_is_hidden: true,
    yax_units: '%',
    target: '#completion',
    x_accessor: 'years',
    y_accessor: 'SE.PRM.CMPT.ZS',
    x_label: 'années',
    mouseover: function(d, i) {
            // custom format the rollover text, show days
            d3.select('#completion svg .mg-active-datapoint')
                .text(d["SE.PRM.CMPT.ZS"] !=null ? '- y: '+ d.years  + ', taux: ' + d["SE.PRM.CMPT.ZS"]: '- y: '+ d.years + ', N/A');
        }
};


var remuneration = {
    title: "Rémunération des enseignants",
    description: "Rémunération de tout le personnel dans l'éducation par année (en fonction du pourcentage des dépenses dans les institutions publiques d'éducation) - Données compilées par la banque mondiale",
    full_width: true,
    height: 350,
    bottom: 65,
    left: 50,
    right: 30,
    min_x: 1999,
    area: false,
    yax_units: '%',
    xax_count: 4,
    missing_is_hidden: true,
    target: '#ensremu',
    x_accessor: 'years',
    y_accessor: 'SE.XPD.MTOT.ZS',
    x_label: 'années',
    mouseover: function(d, i) {
             d3.select('#ensremu svg .mg-active-datapoint')
                .text(d['SE.XPD.MTOT.ZS'] !=null ? '-- y: '+ d.years  + ' ( ' + d['SE.XPD.MTOT.ZS'] +' % )': '-- y: '+ d.years + ' ( N/A )');
        }
};


var westaffcomp = {
    title: "Taux d'achèvement du primaire en 2014",
    description: "Comparaison du taux d'achèvement de l'école primaire en Afrique de l'ouest (2014)",
    full_width: true,
    top: 90,
    height: 350,
    target: '#completion',
    x_accessor: 'years',
    chart_type: 'bar',
    y_accessor: 'country',
    x_accessor: 'value',
        mouseover: function(d, i) {
            // custom format the rollover text, show days
            d3.select('#completion svg .mg-active-datapoint')
                .html(d.country_name + ' ( '+ d.value +' % )');
        }

};

var chooser = {'benin': primary_completion, 'westafr':westaffcomp};


d3.json('/assets/data/total.json', function(data) {
    
    globals.data = data; // MG.convert.date(data, 'years');

    split_by_params.data = data;
    MG.data_graphic(split_by_params);

    primary_completion.data = data;
    MG.data_graphic(primary_completion);

    remuneration.data = data;
    MG.data_graphic(remuneration);


});


d3.json('/assets/data/westafr.json', function(data) {

    westaffcomp.data = data;

});


d3.json('/assets/data/predict.json', function(data) {

    resstagne.data = data;
    MG.data_graphic(resstagne);

});


$('.split-by-controls button').click(function() {
    var new_y_accessor = $(this).data('y_accessor');
    if (new_y_accessor == 'adm'){
        split_by_params.yax_units = '%';
        split_by_params.description = "% d'admis au baccalauréat depuis 1969";

    }else{
        split_by_params.yax_units = '';
        split_by_params.description = "Nombre d'inscrits au baccalauréat depuis 1969";
    }
    split_by_params.y_accessor = new_y_accessor;

    // change button state
    $(this).addClass('active').siblings().removeClass('active');

    // update data
    delete split_by_params.xax_format;
    split_by_params.title = $(this).data('title') + " au bac";
    MG.data_graphic(split_by_params);
});

$('.compselector button').click(function() {
    var show = $(this).data('chooser');
    console.log(show);
    // change button state
    $(this).addClass('active').siblings().removeClass('active');
    MG.data_graphic(chooser[show]);
});

function modify_time_period(data, past_n_years) {
    if (past_n_years !== '') {

        return MG.clone(data).slice(past_n_years * -1);
    }
    return data;
};

function set_marker(graph){
    d3.selectAll('#bacperyear .mg-marker-text')
        .attr('y', 170)
        .style('fill', 'red');
    d3.selectAll('#bacperyear .mg-markers line')
        .attr('y1', 180)
        .attr('y2', 250);

};