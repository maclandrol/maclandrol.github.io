---
layout: page
title: Blog posts
---

<div id="blogContainer" data-page="{{ currentPage }}" data-totalPages="{{ paginator.total_pages }}">
<ul>
{% for post in site.posts limit:site.paginate %} 
	
	<li>{{ post.date | date_to_string }} &raquo; <a href="{{ post.url }}">{{ post.title }}</a>
	</li>
{% endfor %}
	</ul>

</div>

<div class="infinite-spinner"></div>