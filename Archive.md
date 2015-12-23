---
layout: page
title: Blog posts
---

<ul class="post-list">
{% for post in site.posts limit:site.paginate %} 
	<li>{{ post.date | date_to_string }} &raquo; <a href="{{ post.url }}">{{ post.title }}</a>
	</li>
{% endfor %}
</ul>

<div class="infinite-spinner"></div>