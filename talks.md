---
layout: page
title: Talks
---

<ul class="sli-list">
{% for pres in site.data.slides %} 
	<li>{{ pres.date | date_to_string }} : <a href="{{ pres.link }}">{{ pres.title }} - {{ pres.loc }}</a>
	</li>
{% endfor %}
</ul>
