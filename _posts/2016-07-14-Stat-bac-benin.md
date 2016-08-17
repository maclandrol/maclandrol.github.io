---
layout: post
date: 2016-07-14
title: "Sur l'éducation au Bénin, (1/2)."
tags: [stats, statistics, fr]
comments: true
style: ["https://cdnjs.cloudflare.com/ajax/libs/metrics-graphics/2.10.0/metricsgraphics.min.css", "public/css/education.css"]
---


Bénin, __autrefois__ _quartier latin_ de l'Afrique, titre dont certains se pavanent encore fièrement... à tort...

Eh oui, elle est bien loin cette époque. Et pour preuve, les résultats du bac n'ont plus jamais excédés les 50% depuis __1973__, chutant parfois aussi bas que 8.13%. Alors que le nombre d'étudiants particulièrement brillants aux examens scolaires ne cesse de croître ces dernières années, le taux global de réussite aux examens suit, bien malheureusement, la tendance inverse.

<!--more-->

<br>
<div id="bacperyear"> </div>
<div class='button-grp split-by-controls'>
	<button type="button" class="button active" data-y_accessor="adm" data-title="Pourcentage d'admis">Admis</button>
    <button type="button" class="button" data-y_accessor="cand" data-title="Nombre d'inscrits">Inscrits</button>
</div>

39% au Certificat d'Étude Primaire 2016, 16% au Brevet d'Étude du Premier Cycle de la même année. L'on se pose bien des questions.

<div class="midtitle">"Qu'est ce qui n'a pas marché ?"</div>

Le système éducatif du Bénin ne fonctionne pas ! Voilà ce que vous chantera, chaque année, tout béninois déçu par nos piètres résultats aux examens nationaux. Je trouve que cette \-bien vague\- réponse ne fait que mettre en évidence notre ignorance sur les causes potentielles des échecs observés. Est-ce médiocrité des étudiants ? L'incompétence des enseignants ? La faute de l'État ? La réponse doit bien se trouver au milieu de tout ça, n'est ce pas ?



<div class="midtitle">"Faut-il accuser la médiocrité des étudiants"</div>

Vous auriez du mal à me faire gober que les malheureux 51% au CEP sont __tous__ des cancres. Est-ce si difficile que ça d'apprendre à lire et à compter ?  Oui je vous vois venir, donc je reprends : Est-ce si difficile que ça "d'apprendre à lire et à compter en français" ?

Il est évident que le processus d'apprentissage est grandement influencé par la langue d'enseignement. À ce propos, Amadou-Mahtar M'Bow, ex-directeur de l'UNESCO [s'exprimait en ces termes](http://unesdoc.unesco.org/images/0008/000829/082991fo.pdf) :


    ... la langue exprime une certaine vision du monde, est le véhicule
     essentiel d'une culture. C'est elle qui structure les principaux
      processus psychologiques, cognitifs et sociaux par lesquels 
      l'individu s'intégre à la collectivité et se reconnait en elle.
      __L'expérience a prouvé que tous les efforts d'extension de la 
      scolarisation et de l'alphabétisation n'ont qu'un impact limité 
      aussi longtemps qu'ils ne s'intègrent pas dans les réalités de la
       vie quotidienne__, dont la langue constitue un aspect essentiel. 


Voilà qui est dit. Je pense néanmoins, (vous êtes invités à partager mon point de vue) que la langue n'est pas la barrière principale. 

Difficile d'amener un enfant à étudier quand le suivi parental à la maison n'y est pas. Lorsque les parents ne réalisent pas l'importance de l'éducation de leurs enfants, comment voudriez vous que ces derniers soient assidus en classe ? Comment voudriez vous que l'enfant révise ses leçons alors qu'aussitôt rentré à la maison, il doit contribuer à ramener de l'argent parce que sa famille peine à survivre ?

Bien sûr qu'il y a des exceptions, des E-X-C-E-P-T-I-O-N-S. 


L'UNESCO rapportait en 

Devons nous faire confiance à ces stats ? hmmm je vous laisse en débattre...


<div class="midtitle">"Le tier, est ce la limite ?"</div>


<div id="time_period"> </div>
<div class="button-grp">
    <button type="button"  class="button active" data-time_period="">Total</button>
    <button type="button" class="button" data-time_period="10">Décennie passée</button>
</div>

<div class="midtitle"> Le tier, est ce la limite ? </div>

Le système éducatif du Bénin ne fonctionne pas ! Voilà ce que vous chantera, chaque année, tout béninois déçu par la médiocrité des résultats aux examens. Je trouve que cette \-bien vague\- réponse ne fait que mettre en évidence notre ignorance sur les causes potentielles des échecs observés. Est-ce la médiocrité des étudiants ? L'incompétence des enseignants ? La faute de l'État ? La réponse doit bien se trouver au milieu de tout ça. 
% d'alphabete de 2014 au benin des adultes (15 et plus)


"Percentage of teachers in secondary education who are trained, both sexes (%)" "SE.SEC.TCAQ.ZS" "8.7031" en 2014 

"Teachers in secondary education, both sexes (number)"
SE.SEC.TCHR
90439 profs 


<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/d3/4.2.2/d3.min.js"></script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/metrics-graphics/2.10.0/metricsgraphics.js"></script>


<script>

var globals = {};

var split_by_params = {
    title: "Admission au bac",
    description: "% d'admis au bac depuis 1969",
	full_width: true,
    height: 300,
    bottom: 65,
    left: 110,
    right: 40,
    xax_count: 4,
    missing_is_hidden: true,
    missing_text : 'Inscription|Admission au bac depuis 1969',
    target: '#bacperyear',
    x_accessor: 'years',
    show_tooltips: true,
    show_year_markers: true,
    y_accessor: 'adm',
    x_label: 'années',
    y_label: '% admis',
    mouseover: function(d, i) {
            // custom format the rollover text, show days
            d3.select('#bacperyear svg .mg-active-datapoint')
                .text(d.years!=1989 ? 'y: '+ d.years + ' | ' + d.cand +  ' ( ' + d.adm + '% )': '- Année blanche -');
        }


};

var modify_time_period_params = {
    title: "Nombre d'inscrits",
    description: "Nombre d'inscrits par année au baccalauréat béninois",
    full_width: true,
    height: 400,
    right: 40,
    show_secondary_x_label: false,
    xax_count: 4,
    target: '#time_period',
    x_accessor: 'years',
    y_accessor: 'cand'
}

d3.json('http://127.0.0.1:4000/total.json', function(data) {
    
    globals.data = data; // MG.convert.date(data, 'years');

    split_by_params.data = data;

    MG.data_graphic(split_by_params);

    modify_time_period_params.data = data;
    MG.data_graphic(modify_time_period_params);

})

$('.split-by-controls button').click(function() {
    var new_y_accessor = $(this).data('y_accessor');
    if (new_y_accessor == 'adm'){
        split_by_params.y_label  = '% admis';
    }else{
        split_by_params.y_label = 'inscrits';
    }
    split_by_params.y_accessor = new_y_accessor;

    // change button state
    $(this).addClass('active').siblings().removeClass('active');

    // update data
    delete split_by_params.xax_format;
    split_by_params.title = $(this).data('title') + " au bac depuis 1969";
    MG.data_graphic(split_by_params);
});

$('.modify-time-period-controls button').click(function() {
    var past_n_years = $(this).data('time_period');
    var data = modify_time_period(globals.data, past_n_years);

    // change button state
    $(this).addClass('active').siblings().removeClass('active');

    delete modify_time_period_params.xax_format;
    modify_time_period_params.data = data;
    MG.data_graphic(modify_time_period_params);
});

function modify_time_period(data, past_n_years) {
    if (past_n_years !== '') {

        return MG.clone(data).slice(past_n_years * -1);
    }

    return data;
}

function set_marker(graph){
        d3.selectAll('#bacperyear .mg-marker-text')
            .attr('y', 170)
            .style('fill', 'red');
        d3.selectAll('#bacperyear .mg-markers line')
            .attr('y1', 180)
            .attr('y2', 250);

    }
</script>


Current education expenditure, total (% of total expenditure in public institutions) "SE.XPD.CTOT.ZS"
114

"All education staff compensation, total" SE.XPD.MTOT.ZS

118d[]

,"Primary completion rate, total (% of relevant age group)" "SE.PRM.CMPT.ZS"
23

"Expenditure on education as % of total government expenditure (%)" SE.XPD.TOTL.GB.ZS

125

"Government expenditure on education, total (% of GDP)" SE.XPD.TOTL.GD.ZS

